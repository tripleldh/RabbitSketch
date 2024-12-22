#ifndef __PYBIND_H__
#define __PYBIND_H__

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/array.h>
#include <memory>
namespace py = pybind11;



namespace pybind11 { namespace detail {
    template <> struct type_caster<std::unique_ptr<int[]>> {
    public:
        PYBIND11_TYPE_CASTER(std::unique_ptr<int[]>, _("unique_ptr<int[]>"));

        bool load(py::handle src, bool) {
            if (src.is_none()) return false;
            py::array arr = py::cast<py::array>(src);
            auto size = arr.size();
            value = std::make_unique<int[]>(size);
            std::memcpy(value.get(), arr.data(), size * sizeof(int));
            return true;
        }

        static py::handle cast(const std::unique_ptr<int[]>& src, py::return_value_policy, py::handle) {
            auto size = sizeof(src) / sizeof(int);
            return py::array(py::buffer_info(
                src.get(),                        /* Pointer to buffer */
                sizeof(int),                      /* Size of one scalar */
                py::format_descriptor<int>::format(), /* Buffer format */
                1,                                /* Number of dimensions */
                {size},                          /* Shape of the buffer */
                {sizeof(int)}                    /* Strides (in bytes) */
            ));
        }
    };
}}



PYBIND11_MODULE(rabbitsketch, m) {
    m.doc() = "rabbitsketch pybind";
		m.def("read_shuffled_file", &Sketch::read_shuffled_file, "Reads shuffled file and returns a tuple with an array");
  	py::class_<Sketch::MinHash>(m, "MinHash")
  		//.def(py::init<>())
  		//.def(py::init<int>(), py::arg("k"))
  		//.def(py::init<int, int>(), py::arg("k"), py::kwonly(), py::arg("size"))
  		.def(py::init<int, int, uint32_t>(), py::arg("kmer") = 21, py::arg("size")=1000, py::arg("seed")=42)
  		//.def(py::init<int, int, uint32_t>())
  		.def("update", &Sketch::MinHash::update)
  		.def("merge", &Sketch::MinHash::merge)
  		.def("jaccard", &Sketch::MinHash::jaccard)
  		.def("distance", &Sketch::MinHash::distance)
  		.def("get_total_length", &Sketch::MinHash::getTotalLength)
  		.def("print_min_hashes", &Sketch::MinHash::printMinHashes)
  		.def("count", &Sketch::MinHash::count)
  
  		//parameters
  		//.def("setKmerSize", &Sketch::MinHash::setKmerSize)
  		//.def("setAlphabetSize", &Sketch::MinHash::setAlphabetSize)
  		//.def("setPreserveCase", &Sketch::MinHash::setPreserveCase)
  		//.def("setUse64", &Sketch::MinHash::setUse64)
  		//.def("setSeed", &Sketch::MinHash::setSeed)
  		//.def("setSketchSize", &Sketch::MinHash::setSketchSize)
  		//.def("setNoncanonical", &Sketch::MinHash::setNoncanonical)
  		.def("get_kmer_size", &Sketch::MinHash::getKmerSize)
  		.def("get_seed", &Sketch::MinHash::getSeed)
  		.def("get_max_sketch_size", &Sketch::MinHash::getMaxSketchSize)
  		.def("get_sketch_size", &Sketch::MinHash::getSketchSize)
  		.def("is_empty", &Sketch::MinHash::isEmpty)
  		.def("is_reverse_complement", &Sketch::MinHash::isReverseComplement)
  		;
    py::class_<Sketch::OSketch>(m, "OSketch")
        .def(py::init<>())
        .def_readwrite("k", &Sketch::OSketch::k)
        .def_readwrite("l", &Sketch::OSketch::l)
        .def_readwrite("m", &Sketch::OSketch::m)
        .def_readwrite("data", &Sketch::OSketch::data)
        .def_readwrite("rcdata", &Sketch::OSketch::rcdata)
        .def("__eq__", &Sketch::OSketch::operator==);

    py::class_<Sketch::OrderMinHash>(m, "OrderMinHash")
					//.def(py::init<const std::string&>(), py::arg("seqNew"))
    		//.def(py::init<char*>(), py::arg("seqNew"))
				.def(py::init<>())
			  //.def(py::init<char*>())	
				.def("build_sketch", &Sketch::OrderMinHash::buildSketch)
        .def("similarity", &Sketch::OrderMinHash::similarity)
        .def("distance", &Sketch::OrderMinHash::distance)
        .def("set_k", &Sketch::OrderMinHash::setK)
        .def("set_l", &Sketch::OrderMinHash::setL)
        .def("set_m", &Sketch::OrderMinHash::setM)
        .def("set_seed", &Sketch::OrderMinHash::setSeed)
        .def("set_reverse_complement", &Sketch::OrderMinHash::setReverseComplement)
  //      .def("get_k", &Sketch::OrderMinHash::getK)
  //      .def("get_l", &Sketch::OrderMinHash::getL)
  //      .def("get_m", &Sketch::OrderMinHash::getM)
  //      .def("get_seed", &Sketch::OrderMinHash::getSeed)
  //      .def("is_reverse_complement", &Sketch::OrderMinHash::isReverseComplement);
		;

		py::class_<Sketch::HyperLogLog>(m, "HyperLogLog")
						//.def(py::init<>())
						//.def(py::init<int>(), py::arg("k"))
						//.def(py::init<int, int>(), py::arg("k"), py::kwonly(), py::arg("size"))
						.def(py::init<int>(), py::arg("np")=10)
						//.def(py::init<int, int, uint32_t>())
						.def("jaccard", &Sketch::HyperLogLog::jaccard_index)
						//.def("jaccard", py::overload_cast<HyperLogLog&>(&Sketch::HyperLogLog::jaccard_index))
						.def("update", &Sketch::HyperLogLog::update)
						.def("merge", &Sketch::HyperLogLog::merge)
						.def("distance", &Sketch::HyperLogLog::distance)
						;
    py::class_<Sketch::kssd_parameter_t>(m, "kssd_parameter_t")
        .def(py::init<int, int, int, int*>())
        .def_readwrite("half_k", &Sketch::kssd_parameter_t::half_k)
        .def_readwrite("half_subk", &Sketch::kssd_parameter_t::half_subk)
        .def_readwrite("drlevel", &Sketch::kssd_parameter_t::drlevel)
        .def_readwrite("shuffled_dim", &Sketch::kssd_parameter_t::shuffled_dim)
        .def_readwrite("hash_size", &Sketch::kssd_parameter_t::hashSize);
	
    py::class_<Sketch::Kssd>(m, "Kssd")
        .def(py::init<Sketch::kssd_parameter_t>())
        .def("read_shuffled_file", &Sketch::read_shuffled_file, "Reads shuffled file and returns a tuple with an array")
				.def("update", &Sketch::Kssd::update)
        .def("jaccard", &Sketch::Kssd::jaccard)
        .def("distance", &Sketch::Kssd::distance)
        .def("store_hashes", &Sketch::Kssd::storeHashes)
        .def("get_halfk", &Sketch::Kssd::get_half_k)
        .def("get_half_subk", &Sketch::Kssd::get_half_subk)
        .def("get_drlevel", &Sketch::Kssd::get_drlevel);
				
}

#endif // __PYBIND_H__

