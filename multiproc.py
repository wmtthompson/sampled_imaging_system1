from Queue import *
from threading import Thread, Lock

queue = Queue()

threads = 4

def processor():
     if queue.empty() == True:
          print "the Queue is empty!"
          sys.exit(1)
          try:
               job = queue.get()
               print "I'm operating on a job"
               queue.task_done()
          except:
               print "Failed to operate on job"

jobs = ["digging", "shit", "out", "of", "a","hole"]
xx = range(1,5)
for job in jobs:
     queue.put(job)


for i in range(threads):
     th = Thread(target = processor)
     th.setDaemon(True)
     th.start()
queue.join()
