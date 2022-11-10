import random
import math


def get_high_score(f, t):
    if f == t:
        return f
    n = int(math.log2(t))
    #print("max to test 2**", n)
    while n > 0:
        x = int(2**n)
        #print("testing", x, f, t)
        if f - f % x + x <= t:
            #print("hit")
            return f - f % x + x - 1
        n -= 1
    #print("no hit")
    return f


def score(p):
    for x in range(10, 0, -1):
        if ( p + 1 ) % ( 2**x ) == 0:
            return x
    return 0

def find_high_score(f, t):
    s = 0
    x = f
    for i in range(f, t):
        if score(i) > s:
            s = score(i)
            x = i
    return x

PRINT = True

def test(f, t):
    if PRINT:
        for i in range(f, t):
            print(i, end="\t")
        print()
        for i in range(f, t):
            print(score(i), end="\t")
        print()

    a = find_high_score(f, t)
    b = get_high_score(f, t)

    if PRINT:
        print("find_high_score", a, "get_high_score", b)
    if a != b:
        return False
    return True

test(0, 13)
exit()

b = True
while b:
    f = random.randrange(100)
    t = random.randrange(100) + f
    b = test(f, t)