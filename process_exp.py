from multiprocessing import Process
import os

import numpy as np


def info(title):
    print(title)
    print('module name:', __name__)
    print('parent process:', os.getppid())
    print('process id:', os.getpid())

def f(name):
    info('function f')
    print('hello', name)

zzz = list()
def g(x):
    y = x**2
    info("function g")
    print y
    zzz.append(y)

xx = range(1,5)

if __name__ == '__main__':
    info('main line')
    #p = Process(target=f, args=('bob',))
    #p1 = Process(target=f, args=('marcus',))
    for x in xx:
        p = Process(target = g, args=(x,))
        p.start()
        p.join()


    for zz in zzz:
        print "The list y contains"
        print zz

        
''' p3 = Process(target=g, args=(xx[1],))
    p.start()
    p1.start()
    p3.start()
    p.join()
    p1.join()
    p3.join()'''
