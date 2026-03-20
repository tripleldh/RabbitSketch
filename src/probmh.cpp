/**
 * ProbMinHash4 – optimized implementation.
 *
 * Optimizations vs. baseline:
 *  1. Ziggurat exponential distribution  (~5x faster than -log(u))
 *  2. SIMD batched murmur3_fmix  (AVX-512: 8-lane; scalar fallback otherwise)
 *  3. Early return in addHashFromRng BEFORE perm_.reset()
 *  4. Cached class members in hot loop locals
 *  5. Global-max pre-filter: cur_max loaded once per batch (register), refreshed
 *     only when addHashFromRng() is actually called (tighter threshold, zero cost
 *     on pruned lanes)
 *  6. SIMD Jaccard / merge  (AVX-512 / AVX2 / scalar)
 *  7. PermStream version overflow protection
 *  8. No strlen() in update(): caller passes length directly
 *  9. TED parameters packed into TedParam[] (4 arrays → 1, better cache locality)
 * 10. tracker_ leaves ARE the registers – regs_ eliminated (no duplicate storage)
 * 11. O(m) build_from_leaves() for copy / merge (was O(m log m))
 * 12. bool[8] lane_valid array; all-invalid batch skipped via fast OR check
 * 13. Optional single-fmix fast path: compile with -DPMH_FAST_HASH
 *
 * PMH_FAST_HASH notes:
 *   Default (off): two rounds of murmur3_fmix per k-mer (maximum avalanche).
 *   When defined:  one round per k-mer (the SIMD-computed hash is reused as
 *   the RNG seed directly).  Saves ~1 ns/kmer.  Statistical tests (Jaccard
 *   bias, variance, self-similarity) should be verified on real data before
 *   enabling in production.
 */

#include "probmh.h"
#include "hash_int.h"

#include <immintrin.h>
#include <cstring>
#include <cstdio>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>

using namespace Sketch;

// ═══════════════════════════════════════════════════════════════════════════
// RNG utilities
// ═══════════════════════════════════════════════════════════════════════════

static inline uint64_t wymix64(uint64_t a, uint64_t b) {
    __uint128_t r = (__uint128_t)a * b;
    return (uint64_t)(r >> 64) ^ (uint64_t)r;
}

static inline uint64_t wyrand_next(uint64_t& state) {
    state += 0xa0761d6478bd642fULL;
    return wymix64(state, state ^ 0xe7037ed1a0b428dbULL);
}

static inline double wy_uniform(uint64_t& state) {
    return (wyrand_next(state) >> 11) * (1.0 / (1ULL << 53));
}

static inline uint32_t wy_uniform_int(uint32_t n, uint64_t& state) {
    uint64_t x = wyrand_next(state) >> 32;
    uint64_t m = x * static_cast<uint64_t>(n);
    uint32_t l = static_cast<uint32_t>(m);
    if (l < n) {
        uint32_t t = static_cast<uint32_t>(-static_cast<int32_t>(n) % n);
        while (l < t) {
            x = wyrand_next(state) >> 32;
            m = x * static_cast<uint64_t>(n);
            l = static_cast<uint32_t>(m);
        }
    }
    return static_cast<uint32_t>(m >> 32);
}

// ═══════════════════════════════════════════════════════════════════════════
// Ziggurat exponential distribution  (Opt #1)
// Tables from Boost / Ertl's probminhash reference code.
// Fast path (~97%): 1 RNG call + 1 multiply + 1 compare. No log().
// ═══════════════════════════════════════════════════════════════════════════
namespace {

static const double zig_table_x[257] = {
    8.6971174701310497140, 7.6971174701310497140, 6.9410336293772123602,
    6.4783784938325698538, 6.1441646657724730491, 5.8821443157953997963,
    5.6664101674540337371, 5.4828906275260628694, 5.3230905057543986131,
    5.1814872813015010392, 5.0542884899813047117, 4.9387770859012514838,
    4.8329397410251125881, 4.7352429966017412526, 4.6444918854200854873,
    4.5597370617073515513, 4.4802117465284221949, 4.4052876934735729805,
    4.3344436803172730116, 4.2672424802773661873, 4.2033137137351843802,
    4.1423408656640511251, 4.0840513104082974638, 4.0282085446479365106,
    3.9746060666737884793, 3.9230625001354895926, 3.8734176703995089983,
    3.8255294185223367372, 3.7792709924116678992, 3.7345288940397975350,
    3.6912010902374189454, 3.6491955157608538478, 3.6084288131289096339,
    3.5688252656483374051, 3.5303158891293438633, 3.4928376547740601814,
    3.4563328211327607625, 3.4207483572511205323, 3.3860354424603017887,
    3.3521490309001100106, 3.3190474709707487166, 3.2866921715990692095,
    3.2550473085704501813, 3.2240795652862645207, 3.1937579032122407483,
    3.1640533580259734580, 3.1349388580844407393, 3.1063890623398246660,
    3.0783802152540905188, 3.0508900166154554479, 3.0238975044556767713,
    2.9973829495161306949, 2.9713277599210896472, 2.9457143948950456386,
    2.9205262865127406647, 2.8957477686001416838, 2.8713640120155362592,
    2.8473609656351888266, 2.8237253024500354905, 2.8004443702507381944,
    2.7775061464397572041, 2.7548991965623453650, 2.7326126361947007411,
    2.7106360958679293686, 2.6889596887418041593, 2.6675739807732670816,
    2.6464699631518093905, 2.6256390267977886123, 2.6050729387408355373,
    2.5847638202141406911, 2.5647041263169053687, 2.5448866271118700928,
    2.5253043900378279427, 2.5059507635285939648, 2.4868193617402096807,
    2.4679040502973649846, 2.4491989329782498908, 2.4306983392644199088,
    2.4123968126888708336, 2.3942890999214583288, 2.3763701405361408194,
    2.3586350574093374601, 2.3410791477030346875, 2.3236978743901964559,
    2.3064868582835798692, 2.2894418705322694265, 2.2725588255531546952,
    2.2558337743672190441, 2.2392628983129087111, 2.2228425031110364013,
    2.2065690132576635755, 2.1904389667232199235, 2.1744490099377744673,
    2.1585958930438856781, 2.1428764653998416425, 2.1272876713173679737,
    2.1118265460190418108, 2.0964902118017147637, 2.0812758743932248696,
    2.0661808194905755036, 2.0512024094685848641, 2.0363380802487695916,
    2.0215853383189260770, 2.0069417578945183144, 1.9924049782135764992,
    1.9779727009573602295, 1.9636426877895480401, 1.9494127580071845659,
    1.9352807862970511135, 1.9212447005915276767, 1.9073024800183871196,
    1.8934521529393077332, 1.8796917950722108462, 1.8660195276928275962,
    1.8524335159111751661, 1.8389319670188793980, 1.8255131289035192212,
    1.8121752885263901413, 1.7989167704602903934, 1.7857359354841254047,
    1.7726311792313049959, 1.7596009308890742369, 1.7466436519460739352,
    1.7337578349855711926, 1.7209420025219350428, 1.7081947058780575683,
    1.6955145241015377061, 1.6829000629175537544, 1.6703499537164519163,
    1.6578628525741725325, 1.6454374393037234057, 1.6330724165359912048,
    1.6207665088282577216, 1.6085184617988580769, 1.5963270412864831349,
    1.5841910325326886695, 1.5721092393862294810, 1.5600804835278879161,
    1.5481036037145133070, 1.5361774550410318943, 1.5243009082192260050,
    1.5124728488721167573, 1.5006921768428164936, 1.4889578055167456003,
    1.4772686611561334579, 1.4656236822457450411, 1.4540218188487932264,
    1.4424620319720121876, 1.4309432929388794104, 1.4194645827699828254,
    1.4080248915695353509, 1.3966232179170417110, 1.3852585682631217189,
    1.3739299563284902176, 1.3626364025050864742, 1.3513769332583349176,
    1.3401505805295045843, 1.3289563811371163220, 1.3177933761763245480,
    1.3066606104151739482, 1.2955571316866007210, 1.2844819902750125450,
    1.2734342382962410994, 1.2624129290696153434, 1.2514171164808525098,
    1.2404458543344064544, 1.2294981956938491599, 1.2185731922087903071,
    1.2076698934267612830, 1.1967873460884031665, 1.1859245934042023557,
    1.1750806743109117687, 1.1642546227056790397, 1.1534454666557748056,
    1.1426522275816728928, 1.1318739194110786733, 1.1211095477013306083,
    1.1103581087274114281, 1.0996185885325976575, 1.0888899619385472598,
    1.0781711915113727024, 1.0674612264799681530, 1.0567590016025518414,
    1.0460634359770445503, 1.0353734317905289496, 1.0246878730026178052,
    1.0140056239570971074, 1.0033255279156973717, 0.99264640550727647009,
    0.98196705308506317914, 0.97128624098390397896, 0.96060271166866709917,
    0.94991517776407659940, 0.93922231995526297952, 0.92852278474721113999,
    0.91781518207004493915, 0.90709808271569100600, 0.89637001558989069006,
    0.88562946476175228052, 0.87487486629102585352, 0.86410460481100519511,
    0.85331700984237406386, 0.84251035181036928333, 0.83168283773427388393,
    0.82083260655441252290, 0.80995772405741906620, 0.79905617735548788109,
    0.78812586886949324977, 0.77716460975913043936, 0.76617011273543541328,
    0.75513998418198289808, 0.74407171550050873971, 0.73296267358436604916,
    0.72181009030875689912, 0.71061105090965570413, 0.69936248110323266174,
    0.68806113277374858613, 0.67670356802952337911, 0.66528614139267855405,
    0.65380497984766565353, 0.64225596042453703448, 0.63063468493349100113,
    0.61893645139487678178, 0.60715622162030085137, 0.59528858429150359384,
    0.58332771274877027785, 0.57126731653258903915, 0.55910058551154127652,
    0.54682012516331112550, 0.53441788123716615385, 0.52188505159213564105,
    0.50921198244365495319, 0.49638804551867159754, 0.48340149165346224782,
    0.47023927508216945338, 0.45688684093142071279, 0.44332786607355296305,
    0.42954394022541129589, 0.41551416960035700100, 0.40121467889627836229,
    0.38661797794112021568, 0.37169214532991786118, 0.35639976025839443721,
    0.34069648106484979674, 0.32452911701691008547, 0.30783295467493287307,
    0.29052795549123115167, 0.27251318547846547924, 0.25365836338591284433,
    0.23379048305967553619, 0.21267151063096745264, 0.18995868962243277774,
    0.16512762256418831796, 0.13730498094001380420, 0.10483850756582017915,
    0.063852163815003480173, 0
};

static const double zig_table_y[257] = {
    0, 0.00045413435384149675545, 0.00096726928232717452884,
    0.0015362997803015723824, 0.0021459677437189061793,
    0.0027887987935740759640, 0.0034602647778369039855,
    0.0041572951208337952532, 0.0048776559835423925804,
    0.0056196422072054831710, 0.0063819059373191794422,
    0.0071633531836349841425, 0.0079630774380170392396,
    0.0087803149858089752347, 0.0096144136425022094101,
    0.010464810181029979488, 0.011331013597834597488,
    0.012212592426255380661, 0.013109164931254991070,
    0.014020391403181937334, 0.014945968011691148079,
    0.015885621839973162490, 0.016839106826039946359,
    0.017806200410911360563, 0.018786700744696029497,
    0.019780424338009741737, 0.020787204072578117603,
    0.021806887504283582125, 0.022839335406385238829,
    0.023884420511558170348, 0.024942026419731782971,
    0.026012046645134218076, 0.027094383780955798424,
    0.028188948763978634421, 0.029295660224637394015,
    0.030414443910466605492, 0.031545232172893605499,
    0.032687963508959533317, 0.033842582150874329031,
    0.035009037697397411067, 0.036187284781931419754,
    0.037377282772959360128, 0.038578995503074859626,
    0.039792391023374122670, 0.041017441380414820816,
    0.042254122413316231413, 0.043502413568888183301,
    0.044762297732943280694, 0.046033761076175166762,
    0.047316792913181548703, 0.048611385573379494401,
    0.049917534282706374944, 0.051235237055126279830,
    0.052564494593071689595, 0.053905310196046085104,
    0.055257689676697038322, 0.056621641283742874438,
    0.057997175631200659098, 0.059384305633420264487,
    0.060783046445479636051, 0.062193415408540996150,
    0.063615431999807331076, 0.065049117786753755036,
    0.066494496385339779043, 0.067951593421936607770,
    0.069420436498728751675, 0.070901055162371828426,
    0.072393480875708743023, 0.073897746992364746308,
    0.075413888734058408453, 0.076941943170480510100,
    0.078481949201606426042, 0.080033947542319910023,
    0.081597980709237420930, 0.083174093009632380354,
    0.084762330532368125386, 0.086362741140756912277,
    0.087975374467270219300, 0.089600281910032864534,
    0.091237516631040162057, 0.092887133556043546523,
    0.094549189376055853718, 0.096223742550432800103,
    0.097910853311492199618, 0.099610583670637128826,
    0.10132299742595363588, 0.10304816017125771553,
    0.10478613930657016928, 0.10653700405000166218,
    0.10830082545103379867, 0.11007767640518539026,
    0.11186763167005629731, 0.11367076788274431301,
    0.11548716357863353664, 0.11731689921155557057,
    0.11916005717532768467, 0.12101672182667483729,
    0.12288697950954513498, 0.12477091858083096578,
    0.12666862943751066518, 0.12858020454522817870,
    0.13050573846833078225, 0.13244532790138752023,
    0.13439907170221363078, 0.13636707092642885841,
    0.13834942886358021406, 0.14034625107486244210,
    0.14235764543247220043, 0.14438372216063476473,
    0.14642459387834493787, 0.14848037564386679222,
    0.15055118500103990354, 0.15263714202744286154,
    0.15473836938446807312, 0.15685499236936522013,
    0.15898713896931420572, 0.16113493991759203183,
    0.16329852875190180795, 0.16547804187493600915,
    0.16767361861725019322, 0.16988540130252766513,
    0.17211353531532005700, 0.17435816917135348788,
    0.17661945459049489581, 0.17889754657247831241,
    0.18119260347549629488, 0.18350478709776746150,
    0.18583426276219711495, 0.18818119940425430485,
    0.19054576966319540013, 0.19292814997677133873,
    0.19532852067956322315, 0.19774706610509886464,
    0.20018397469191127727, 0.20263943909370901930,
    0.20511365629383770880, 0.20760682772422204205,
    0.21011915938898825914, 0.21265086199297827522,
    0.21520215107537867786, 0.21777324714870053264,
    0.22036437584335949720, 0.22297576805812018050,
    0.22560766011668406495, 0.22826029393071670664,
    0.23093391716962742173, 0.23362878343743333945,
    0.23634515245705964715, 0.23908329026244917002,
    0.24184346939887722761, 0.24462596913189210901,
    0.24743107566532763894, 0.25025908236886230967,
    0.25311029001562948171, 0.25598500703041538015,
    0.25888354974901621678, 0.26180624268936295243,
    0.26475341883506220209, 0.26772541993204481808,
    0.27072259679906003167, 0.27374530965280298302,
    0.27679392844851734458, 0.27986883323697289920,
    0.28297041453878076010, 0.28609907373707684673,
    0.28925522348967773308, 0.29243928816189258772,
    0.29565170428126120948, 0.29889292101558177099,
    0.30216340067569352897, 0.30546361924459023541,
    0.30879406693456016794, 0.31215524877417956945,
    0.31554768522712893632, 0.31897191284495723773,
    0.32242848495608914289, 0.32591797239355619822,
    0.32944096426413633091, 0.33299806876180896713,
    0.33658991402867758144, 0.34021714906678004560,
    0.34388044470450243010, 0.34758049462163698567,
    0.35131801643748334681, 0.35509375286678745925,
    0.35890847294874976196, 0.36276297335481777335,
    0.36665807978151414890, 0.37059464843514599421,
    0.37457356761590215193, 0.37859575940958081092,
    0.38266218149600982112, 0.38677382908413768115,
    0.39093173698479710717, 0.39513698183329015336,
    0.39939068447523107877, 0.40369401253053026739,
    0.40804818315203238238, 0.41245446599716116772,
    0.41691418643300289465, 0.42142872899761659635,
    0.42599954114303435739, 0.43062813728845883923,
    0.43531610321563659758, 0.44006510084235387501,
    0.44487687341454851593, 0.44975325116275498919,
    0.45469615747461548049, 0.45970761564213768669,
    0.46478975625042618067, 0.46994482528395999841,
    0.47517519303737738299, 0.48048336393045423016,
    0.48587198734188493564, 0.49134386959403255500,
    0.49690198724154955294, 0.50254950184134769289,
    0.50828977641064283495, 0.51412639381474855788,
    0.52006317736823356823, 0.52610421398361972602,
    0.53225388026304326945, 0.53851687200286186590,
    0.54489823767243963663, 0.55140341654064131685,
    0.55803828226258748140, 0.56480919291240022434,
    0.57172304866482579008, 0.57878735860284503057,
    0.58601031847726802755, 0.59340090169173341521,
    0.60096896636523224742, 0.60872538207962206507,
    0.61668218091520762326, 0.62485273870366592605,
    0.63325199421436607968, 0.64189671642726607018,
    0.65080583341457104881, 0.66000084107899974178,
    0.66950631673192477684, 0.67935057226476538741,
    0.68956649611707798890, 0.70019265508278816709,
    0.71127476080507597882, 0.72286765959357200702,
    0.73503809243142351530, 0.74786862198519510742,
    0.76146338884989624862, 0.77595685204011559675,
    0.79152763697249565519, 0.80842165152300838005,
    0.82699329664305033399, 0.84778550062398962096,
    0.87170433238120363669, 0.90046992992574643800,
    0.93814368086217467916, 1
};

} // anon namespace

// Ziggurat Exp(1): fast path needs no log() call (~97% of the time).
static inline double zig_exponential(uint64_t& state) {
    double shift = 0;
    for (;;) {
        uint64_t bits = wyrand_next(state);
        int i = (int)(bits & 0xFF);
        double u = (double)(bits >> 11) * (1.0 / (1ULL << 53));
        double x = u * zig_table_x[i];
        if (x < zig_table_x[i + 1]) return shift + x;
        if (i == 0) { shift += zig_table_x[1]; continue; }
        double y01 = wy_uniform(state);
        double y = zig_table_y[i] + y01 * (zig_table_y[i+1] - zig_table_y[i]);
        double y_above_ubound = (zig_table_x[i] - zig_table_x[i+1]) * y01
                                - (zig_table_x[i] - x);
        double y_above_lbound = y - (zig_table_y[i+1]
                                + (zig_table_x[i+1] - x) * zig_table_y[i+1]);
        if (y_above_ubound < 0 &&
            (y_above_lbound < 0 || y < std::exp(-x)))
            return x + shift;
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Truncated-exponential sampler  (unchanged, fast path dominates)
// ═══════════════════════════════════════════════════════════════════════════

static inline double ted_sample(double c1, double c2, double c3,
                                uint64_t& rng) {
    double x = wy_uniform(rng) * c1;
    if (x < 1.0) return x;
    while (true) {
        x = wy_uniform(rng);
        if (x <= c2) return x;
        double y = wy_uniform(rng) * 0.5;
        if (y > 1.0 - x) { x = 1.0 - x; y = 1.0 - y; }
        if (x <= c3 * (1.0 - y)) return x;
        double c1y = c1 * y;
        if (c1y <= 1.0 - x) return x;
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// DNA encoding  (same LUT as SetSketch)
// ═══════════════════════════════════════════════════════════════════════════

static const uint8_t ENC_LUT[256] = {
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,  0,255,  1,255,255,255,  2,255,255,255,255,255,255,255,255,
    255,255,255,255,  3,255,255,255,255,255,255,255,255,255,255,255,
    255,  0,255,  1,255,255,255,  2,255,255,255,255,255,255,255,255,
    255,255,255,255,  3,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
};
#define PMH_ENC(c)   (ENC_LUT[(uint8_t)(c)])
#define PMH_COMP(e)  ((uint8_t)((e) <= 3 ? 3 - (e) : 255))
#define PMH_VALID(e) ((e) <= 3)

static const double PMH_INV2_53 = 1.0 / (double)(1ULL << 53);


// ═══════════════════════════════════════════════════════════════════════════
// swap helpers
// ═══════════════════════════════════════════════════════════════════════════

void Sketch::swap(ProbMHMaxTracker& a, ProbMHMaxTracker& b) noexcept {
    std::swap(a.m_,       b.m_);
    std::swap(a.lastIdx_, b.lastIdx_);
    std::swap(a.v_,       b.v_);
}

void Sketch::swap(ProbMHPermStream& a, ProbMHPermStream& b) noexcept {
    std::swap(a.m_,       b.m_);
    std::swap(a.idx_,     b.idx_);
    std::swap(a.ver_,     b.ver_);
    std::swap(a.val_,     b.val_);
    std::swap(a.ver_arr_, b.ver_arr_);
}

// ═══════════════════════════════════════════════════════════════════════════
// ProbMHMaxTracker
// ═══════════════════════════════════════════════════════════════════════════

ProbMHMaxTracker::ProbMHMaxTracker(uint32_t m)
    : m_(m), lastIdx_((m << 1) - 2), v_(new double[(m << 1) - 1]) {}

void ProbMHMaxTracker::reset(double inf) {
    std::fill_n(v_.get(), lastIdx_ + 1, inf);
}

bool ProbMHMaxTracker::update(uint32_t idx, double value) {
    assert(idx < m_);
    if (value >= v_[idx]) return false;
    while (true) {
        v_[idx] = value;
        const uint32_t par = m_ + (idx >> 1);
        if (par > lastIdx_) break;
        const uint32_t sib = idx ^ UINT32_C(1);
        const double   sv  = v_[sib];
        if (!(sv < v_[par])) break;
        if (value < sv) value = sv;
        idx = par;
    }
    return true;
}

bool ProbMHMaxTracker::isUpdatePossible(double value) const {
    return value < v_[lastIdx_];
}

double ProbMHMaxTracker::getMax() const {
    return v_[lastIdx_];
}

// O(m) bottom-up rebuild of all internal nodes from pre-set leaves.
// For internal node p, children are at (p-m)*2 and (p-m)*2+1.
// Both children have indices < p, so left-to-right sweep is correct.
void ProbMHMaxTracker::build_from_leaves() {
    for (uint32_t p = m_; p <= lastIdx_; ++p) {
        const uint32_t l = (p - m_) << 1;
        v_[p] = std::max(v_[l], v_[l + 1]);
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// ProbMHPermStream  (Opt #7: overflow-safe version counter)
// ═══════════════════════════════════════════════════════════════════════════

ProbMHPermStream::ProbMHPermStream(uint32_t m)
    : m_(m), idx_(0), ver_(0),
      val_(new uint32_t[m]),
      ver_arr_(new uint32_t[m]) {}

void ProbMHPermStream::reset() {
    idx_ = 0;
    if (ver_ == 0 || ver_ == UINT32_MAX) {
        for (uint32_t i = 0; i < m_; ++i) { val_[i] = i; ver_arr_[i] = 0; }
        ver_ = 1;
    } else {
        ++ver_;
    }
}

uint32_t ProbMHPermStream::next(uint64_t& rng) {
    const uint32_t k = idx_ + wy_uniform_int(m_ - idx_, rng);
    const uint32_t v = ver_;
    const uint32_t result = (ver_arr_[k] != v) ? k : val_[k];
    const uint32_t x     = (ver_arr_[idx_] != v) ? idx_ : val_[idx_];
    val_[k]     = x;
    ver_arr_[k] = v;
    ++idx_;
    return result;
}

// ═══════════════════════════════════════════════════════════════════════════
// ProbMinHash4  –  constructor / copy / assign
// ═══════════════════════════════════════════════════════════════════════════

ProbMinHash4::ProbMinHash4(uint32_t m, int kmer_size, uint64_t seed)
    : m_(m), kmer_size_(kmer_size), seed_(seed),
      ted_params_(new TedParam[m - 1]),
      tracker_(m),
      perm_(m)
{
    assert(m > 1);
    assert(kmer_size >= 1 && kmer_size <= 32);
    tracker_.reset(std::numeric_limits<double>::infinity());

    const double firstBoundary = std::log1p(1.0 / static_cast<double>(m - 1));
    firstBoundaryInv_ = 1.0 / firstBoundary;

    // Index 0: boundary = 1.0, gap unused in delta formula.
    {
        const double rate = firstBoundary;
        ted_params_[0].boundary = 1.0;
        ted_params_[0].gap      = 1.0;
        ted_params_[0].c1 = (rate != 0.0) ? std::expm1(rate) / rate : 1.0;
        ted_params_[0].c2 = (rate != 0.0) ? -std::log1p(std::expm1(-rate) * 0.5) / rate : 0.5;
        ted_params_[0].c3 = (rate != 0.0) ? -std::expm1(-rate) / rate : 1.0;
    }

    double prevBoundary = firstBoundary;
    for (uint32_t i = 1; i < m - 1; ++i) {
        const double b    = std::log1p(static_cast<double>(i + 1) /
                                       static_cast<double>(m - i - 1));
        const double rate = b - prevBoundary;
        const double bNorm = b / firstBoundary;
        ted_params_[i].boundary = bNorm;
        ted_params_[i].gap      = bNorm - ted_params_[i - 1].boundary;
        ted_params_[i].c1 = (rate != 0.0) ? std::expm1(rate) / rate : 1.0;
        ted_params_[i].c2 = (rate != 0.0) ? -std::log1p(std::expm1(-rate) * 0.5) / rate : 0.5;
        ted_params_[i].c3 = (rate != 0.0) ? -std::expm1(-rate) / rate : 1.0;
        prevBoundary = b;
    }
}

ProbMinHash4::ProbMinHash4(const ProbMinHash4& o)
    : m_(o.m_), kmer_size_(o.kmer_size_), seed_(o.seed_),
      ted_params_(new TedParam[o.m_ - 1]),
      firstBoundaryInv_(o.firstBoundaryInv_),
      tracker_(o.m_),
      perm_(o.m_)
{
    std::copy(o.ted_params_.get(), o.ted_params_.get() + m_ - 1,
              ted_params_.get());
    // Copy leaves then rebuild internal nodes in O(m) instead of m × update().
    std::copy(o.tracker_.leaves(), o.tracker_.leaves() + m_,
              tracker_.leaves());
    tracker_.build_from_leaves();
}

ProbMinHash4& ProbMinHash4::operator=(ProbMinHash4 other) {
    std::swap(m_,                other.m_);
    std::swap(kmer_size_,        other.kmer_size_);
    std::swap(seed_,             other.seed_);
    std::swap(ted_params_,       other.ted_params_);
    std::swap(firstBoundaryInv_, other.firstBoundaryInv_);
    Sketch::swap(tracker_,       other.tracker_);
    Sketch::swap(perm_,          other.perm_);
    return *this;
}

// ═══════════════════════════════════════════════════════════════════════════
// addHash  (Opt #3: early return BEFORE perm_.reset())
// ═══════════════════════════════════════════════════════════════════════════

void ProbMinHash4::addHash(uint64_t h) {
    addHashFromRng(mc::murmur3_fmix(h, seed_));
}

void ProbMinHash4::addHashFromRng(uint64_t rng) {
    // Fast path (~97%): first sample u*c1[0] < 1.
    const TedParam& tp0 = ted_params_[0];
    const double u0  = (double)(rng >> 11) * PMH_INV2_53;
    const double hv0 = u0 * tp0.c1;
    double hv;
    if (__builtin_expect(hv0 < 1.0, 1))
        hv = hv0;
    else
        hv = ted_sample(tp0.c1, tp0.c2, tp0.c3, rng);

    if (!tracker_.isUpdatePossible(hv)) return;

    perm_.reset();

    uint32_t i = 1;
    while (tracker_.isUpdatePossible(hv)) {
        tracker_.update(perm_.next(rng), hv);
        if (!tracker_.isUpdatePossible(ted_params_[i - 1].boundary)) break;
        if (i < m_ - 1) {
            const TedParam& tp = ted_params_[i];
            hv = ted_params_[i - 1].boundary +
                 tp.gap * ted_sample(tp.c1, tp.c2, tp.c3, rng);
        } else {
            hv = ted_params_[m_ - 2].boundary +
                 firstBoundaryInv_ * zig_exponential(rng);
            if (tracker_.isUpdatePossible(hv))
                tracker_.update(perm_.next(rng), hv);
            break;
        }
        ++i;
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// update  (Opts #2,4,5,8,12,13)
// ═══════════════════════════════════════════════════════════════════════════

void ProbMinHash4::update(const char* seq, uint64_t length) {
    const int K = kmer_size_;
    if (length < static_cast<uint64_t>(K)) return;

    const uint64_t loc_seed = seed_;

    const uint64_t kmer_mask = (K == 32) ? ~0ULL : ((1ULL << (2 * K)) - 1);
    uint64_t fwd = 0, rev = 0;
    int inv = 0;
    for (int k = 0; k < K; ++k) {
        uint8_t ef = PMH_ENC(seq[k]);
        if (!PMH_VALID(ef)) inv++;
        fwd = (fwd << 2) | (PMH_VALID(ef) ? (ef & 3u) : 0u);
        uint8_t er = PMH_VALID(ef) ? PMH_COMP(ef) : 0u;
        rev = (rev >> 2) | (static_cast<uint64_t>(er & 3u) << (2 * (K - 1)));
    }

    const int      lanes    = 8;
    const uint64_t N_body   = length - static_cast<uint64_t>(K);
    const uint64_t N_batch  = (N_body >= 1) ? (N_body / lanes) * lanes : 0;
    const double   loc_c1_0 = ted_params_[0].c1;

    for (uint64_t i = 0; i < N_batch; i += lanes) {
        uint64_t resv[8];
        bool     lane_valid[8];

        for (int j = 0; j < lanes; ++j) {
            uint64_t pos = i + j;
            lane_valid[j] = (inv == 0);
            resv[j] = lane_valid[j] ? ((fwd <= rev) ? fwd : rev) : 0;
            uint8_t ef_out = PMH_ENC(seq[pos]);
            uint8_t ef_in  = PMH_ENC(seq[pos + K]);
            if (!PMH_VALID(ef_out)) inv--;
            if (!PMH_VALID(ef_in))  inv++;
            fwd = ((fwd << 2) | (PMH_VALID(ef_in) ? (ef_in & 3u) : 0u)) & kmer_mask;
            uint8_t er_in = PMH_VALID(ef_in) ? (PMH_COMP(ef_in) & 3u) : 0u;
            rev = (rev >> 2) | (static_cast<uint64_t>(er_in) << (2 * (K - 1)));
        }

        // Skip SIMD hash if every k-mer in this batch contains N.
        if (!(lane_valid[0]|lane_valid[1]|lane_valid[2]|lane_valid[3]|
              lane_valid[4]|lane_valid[5]|lane_valid[6]|lane_valid[7])) continue;

        // ── Opt #2: 8-lane batched murmur3_fmix (AVX-512) ───────────────
        uint64_t hashvalv[8];
#if defined(__AVX512F__) && defined(__AVX512DQ__)
        {
            __m512i vb = _mm512_loadu_si512((const void*)resv);
            __m512i vs = _mm512_set1_epi64((int64_t)loc_seed);
            __m512i va = _mm512_xor_epi64(vb, vs);
            __m512i vt = _mm512_srli_epi64(va, 33);
            vb = _mm512_xor_epi64(va, vt);
            va = _mm512_mullo_epi64(vb, _mm512_set1_epi64(0xff51afd7ed558ccdLL));
            vt = _mm512_srli_epi64(va, 33);
            vb = _mm512_xor_epi64(va, vt);
            va = _mm512_mullo_epi64(vb, _mm512_set1_epi64(0xc4ceb9fe1a85ec53LL));
            vt = _mm512_srli_epi64(va, 33);
            vb = _mm512_xor_epi64(va, vt);
            _mm512_storeu_si512(hashvalv, vb);
        }
#else
        for (int j = 0; j < lanes; ++j)
            hashvalv[j] = mc::murmur3_fmix(resv[j], loc_seed);
#endif

        // ── Opt #5: cur_max loaded once per batch into a register.
        // Refreshed only when addHashFromRng() is actually invoked, so
        // later lanes in the same batch benefit from the tighter threshold
        // without paying a memory load on every pruned lane.
        double cur_max = tracker_.getMax();
        for (int j = 0; j < lanes; ++j) {
            if (!lane_valid[j]) continue;
#ifdef PMH_FAST_HASH
            uint64_t rng_inner = hashvalv[j];
#else
            uint64_t rng_inner = mc::murmur3_fmix(hashvalv[j], loc_seed);
#endif
            double hv_fast = (double)(rng_inner >> 11) * PMH_INV2_53 * loc_c1_0;
            if (__builtin_expect(hv_fast < 1.0, 1) && hv_fast >= cur_max) continue;
            addHashFromRng(rng_inner);
            cur_max = tracker_.getMax();   // refresh after real update attempt
        }
    }

    // ── Remainder loop (scalar) ──────────────────────────────────────────
    for (uint64_t i = N_batch; i <= N_body; ++i) {
        if (inv == 0) {
            const uint64_t canonical = (fwd <= rev) ? fwd : rev;
#ifdef PMH_FAST_HASH
            uint64_t rng_inner = mc::murmur3_fmix(canonical, loc_seed);
#else
            uint64_t h         = mc::murmur3_fmix(canonical, loc_seed);
            uint64_t rng_inner = mc::murmur3_fmix(h, loc_seed);
#endif
            double hv_fast = (double)(rng_inner >> 11) * PMH_INV2_53 * loc_c1_0;
            if (!(__builtin_expect(hv_fast < 1.0, 1) &&
                  hv_fast >= tracker_.getMax()))
                addHashFromRng(rng_inner);
        }
        if (i < N_body) {
            uint8_t ef_out = PMH_ENC(seq[i]);
            uint8_t ef_in  = PMH_ENC(seq[i + K]);
            if (!PMH_VALID(ef_out)) inv--;
            if (!PMH_VALID(ef_in))  inv++;
            fwd = ((fwd << 2) | (PMH_VALID(ef_in) ? (ef_in & 3u) : 0u)) & kmer_mask;
            uint8_t er_in = PMH_VALID(ef_in) ? (PMH_COMP(ef_in) & 3u) : 0u;
            rev = (rev >> 2) | (static_cast<uint64_t>(er_in) << (2 * (K - 1)));
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// jaccard  (Opt #6: SIMD comparison — AVX-512 / AVX2 / scalar)
// ═══════════════════════════════════════════════════════════════════════════

double ProbMinHash4::jaccard(const ProbMinHash4& other) const {
    assert(m_ == other.m_);
    const double* __restrict__ a = tracker_.leaves();
    const double* __restrict__ b = other.tracker_.leaves();
    const double inf = std::numeric_limits<double>::infinity();
    int count = 0;
    uint32_t k = 0;

#if defined(__AVX512F__)
    {
        __m512d vinf = _mm512_set1_pd(inf);
        for (; k + 8 <= m_; k += 8) {
            __m512d va = _mm512_loadu_pd(a + k);
            __m512d vb = _mm512_loadu_pd(b + k);
            __mmask8 eq     = _mm512_cmp_pd_mask(va, vb, _CMP_EQ_OQ);
            __mmask8 notinf = _mm512_cmp_pd_mask(va, vinf, _CMP_NEQ_UQ);
            count += __builtin_popcount(eq & notinf);
        }
    }
#elif defined(__AVX2__)
    {
        __m256d vinf = _mm256_set1_pd(inf);
        for (; k + 4 <= m_; k += 4) {
            __m256d va  = _mm256_loadu_pd(a + k);
            __m256d vb  = _mm256_loadu_pd(b + k);
            __m256d ceq = _mm256_cmp_pd(va, vb, _CMP_EQ_OQ);
            __m256d cni = _mm256_cmp_pd(va, vinf, _CMP_NEQ_UQ);
            __m256d res = _mm256_and_pd(ceq, cni);
            count += __builtin_popcount(_mm256_movemask_pd(res));
        }
    }
#endif

    for (; k < m_; ++k)
        count += (a[k] == b[k] && a[k] != inf) ? 1 : 0;

    return static_cast<double>(count) / static_cast<double>(m_);
}

// ═══════════════════════════════════════════════════════════════════════════
// merge  –  element-wise min; O(m) tracker rebuild (was O(m log m))
// ═══════════════════════════════════════════════════════════════════════════

ProbMinHash4 ProbMinHash4::merge(const ProbMinHash4& other) const {
    assert(m_ == other.m_ && kmer_size_ == other.kmer_size_);
    ProbMinHash4 ret(m_, kmer_size_, seed_);

    const double* __restrict__ a = tracker_.leaves();
    const double* __restrict__ b = other.tracker_.leaves();
    double*       __restrict__ r = ret.tracker_.leaves();
    uint32_t k = 0;

#if defined(__AVX512F__)
    for (; k + 8 <= m_; k += 8) {
        __m512d vr = _mm512_min_pd(_mm512_loadu_pd(a + k),
                                   _mm512_loadu_pd(b + k));
        _mm512_storeu_pd(r + k, vr);
    }
#elif defined(__AVX2__)
    for (; k + 4 <= m_; k += 4) {
        __m256d vr = _mm256_min_pd(_mm256_loadu_pd(a + k),
                                   _mm256_loadu_pd(b + k));
        _mm256_storeu_pd(r + k, vr);
    }
#endif
    for (; k < m_; ++k)
        r[k] = std::min(a[k], b[k]);

    // Single O(m) sweep to build the tournament-tree internal nodes.
    ret.tracker_.build_from_leaves();
    return ret;
}

// ═══════════════════════════════════════════════════════════════════════════
// print
// ═══════════════════════════════════════════════════════════════════════════

void ProbMinHash4::printSketch() const {
    const double* regs = tracker_.leaves();
    std::fprintf(stdout, "ProbMinHash4 m=%u k=%d regs[0..19]: ", m_, kmer_size_);
    for (uint32_t i = 0; i < m_ && i < 20; ++i)
        std::fprintf(stdout, "%.4g ", regs[i]);
    if (m_ > 20) std::fprintf(stdout, "...");
    std::fprintf(stdout, "\n");
}

#undef PMH_ENC
#undef PMH_COMP
#undef PMH_VALID
