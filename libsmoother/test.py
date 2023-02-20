from threading import Thread
from time import sleep
from libsmoother import test


def print_stuff():
    for _ in range(3):
        print("stuff")
        sleep(1)


def test_2():
    print("sleeping py...")
    sleep(3)
    print("done sleeping py...")


def cm_test():
    thread = Thread(target=test)
    thread.start()
    return thread


def cm_test_2():
    thread = Thread(target=test_2)
    thread.start()
    return thread


print("A")
test()

print("B")
print_stuff()


print("C")
t = cm_test()
print_stuff()
t.join()

print("D")
t = cm_test_2()
print_stuff()
t.join()

print("E")
