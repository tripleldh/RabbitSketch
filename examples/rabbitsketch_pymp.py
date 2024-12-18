import rabbitsketch as sketch
import fastx,time  
import pymp
from multiprocessing import Manager,Queue
import math
import numpy as np
import  dill as pickle
from dill import dumps, loads
file_list = []
t1 = time.time()
manager = Manager()
with open('bacteria.list', "r") as bact1000:
    for genomefile in bact1000:
        genomefile = genomefile.strip()
        file_list.append(genomefile)
print("length of file is:", len(file_list))
def write_to_file(file_name, content):
    with open(file_name, 'a', encoding='utf-8') as file:
        file.write(content)
with pymp.Parallel(128) as p:
    for index in p.xrange(0, len(file_list)):
            kssd = sketch.MinHash()
            for name, seq, qual in fastx.Fastx(file_list[index]):
                kssd.update(seq)
#            simple save sketch case
#            write_to_file(str(p.thread_num)+".txt", kssd.printMinHashes())
t2 = time.time()
print("sketch time is :", t2-t1)
