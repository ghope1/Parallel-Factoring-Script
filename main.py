import random
import time
import math
import fractions
import multiprocessing

from multiprocessing import Process

N = 2

# this function performs the extended euclid algorithm recursively
# this is essentially verbatim the pseudocode from class
def extended_euclid(a, b):
    if b == 0:
        return a, 1, 0
    else:
        old_d, old_x, old_y = extended_euclid(b, a % b)

    d, x, y = old_d, old_y, old_x - math.floor(a//b)*old_y
    return d, x, y


# this function performs a brute force search
def bruteForce(numberToFactor):
    # start by checking initial 2-100
    for i in range(2, min(numberToFactor -1, 100)):
        if numberToFactor % i == 0:
            print(numberToFactor, "has a factor of:", i)
            return i
    # start real search from sqrt(numberToFactor)
    currCheck = math.ceil(pow(numberToFactor, 1/2))

    # if number being checked is even, add one to make it odd
    if currCheck % 2 == 0:
        currCheck += 1

    # keep checking until a factor is found
    while numberToFactor % currCheck != 0:
        # subtract by 2 so you only check odd numbers
        currCheck -= 2
    print(numberToFactor, "has a factor of:",  currCheck)
    return currCheck


# this function performs a search i found online here: https://arxiv.org/pdf/1408.2608.pdf
# i honestly have no clue why it works i just implemented it based on the description in the document
# it works but is pretty slow so i didn't use it over pollard-rho
def pollardStrassen(numberToFactor):
    x0 = min(math.ceil(pow(17*numberToFactor, 1/3)), math.floor(pow(numberToFactor, 1/2)))
    x = x0 + 2
    H = 1

    for i in range(2, x0+1):
        if numberToFactor % i == 0:
            return i
    while x - H <= math.floor(math.pow(numberToFactor, 1/2)):
        cf_convs = continued_fraction_convergents(numberToFactor//(pow(x, 2)))
        maximum_den = 1
        maximum_num = cf_convs[1]
        for tests in cf_convs:
            test = fractions.Fraction(tests).limit_denominator(10000)
            if 4*H > test.denominator > maximum_den:
                maximum_den = test.denominator
                maximum_num = test.numerator

        b = maximum_num
        q = maximum_den
        a = math.floor((q*numberToFactor)/x)

        if math.pow((b * x - a), 2) - 4 * (q * numberToFactor - a * x) * b >= 0 and b != 0:
            h1 = -((b*x-a) + math.pow(math.pow((b*x-a), 2)-4*(q*numberToFactor-a*x)*b, 1/2))/(2*b)
            h2 = -((b*x-a) - math.pow(math.pow((b*x-a), 2)-4*(q*numberToFactor-a*x)*b, 1/2))/(2*b)
            if numberToFactor % (x+h1) == 0:
                if (x+h1) != math.floor(x+h1):
                    print(numberToFactor, "has a factor of:", (x+h1)*fractions.Fraction(x+h1).limit_denominator(1000000).denominator)
                    return (x+h1)*fractions.Fraction(x+h1).limit_denominator(1000000).denominator
                return x + h1
            if numberToFactor % (x + h2) == 0:
                if (x+h2) != math.floor(x+h2):
                    print(numberToFactor, "has a factor of:", (x+h2)*fractions.Fraction(x+h2).limit_denominator(1000000).denominator)
                    return (x+h2)*fractions.Fraction(x+h2).limit_denominator(1000000).denominator
                return x + h2
        x = x+2*H+1
        H = math.floor(math.pow(17*numberToFactor, -1/3)*x)


# this function is required for pollardStrassen to work
# it returns the first 40 continued fraction convergents
# it can be used to produce a fraction that converges to an irrational number
def continued_fraction_convergents(number):
    temp = number

    a_vals = [math.floor(temp)]*40
    cons = [math.floor(temp)]*40

    for i in range(40):
        a_vals[i] = math.floor(temp)
        temp = 1/(temp-a_vals[i])
        if i == 0:
            cons[i] = a_vals[i]

        elif i == 1:
            cons[i] = (1 + a_vals[0]*a_vals[1])/a_vals[1]

        else:
            minus_1 = fractions.Fraction(cons[i - 1]).limit_denominator(1000000)
            minus_2 = fractions.Fraction(cons[i - 2]).limit_denominator(1000000)
            if (a_vals[i]*minus_1.denominator+minus_2.denominator) != 0:
                cons[i] = (a_vals[i]*minus_1.numerator+minus_2.numerator)/(a_vals[i]*minus_1.denominator+minus_2.denominator)

    return cons


# this function performs the pollard-rho search based on the pseudocode from the typed notes from class
def pollardRho(numberToFactor, lowerBound, upperBound):
    i = 1
    currX = random.randint(lowerBound, upperBound)
    y = currX
    k = 2
    d = 1

    # this does not use while(true) like the typed notes do
    while d == 1 or d == numberToFactor:
        i = i + 1
        prevX = currX
        currX = (pow(prevX, 2) - 1) % numberToFactor

        d = extended_euclid(y - currX, numberToFactor)

        d = d[0]

        if d != 1 and d != numberToFactor:
            print(numberToFactor, "has a factor of:", d)
        if i == k:
            y = currX
            k = 2*k
    return d


def main():
    # number of processes running pollard-rho
    numPollardRho = 7

    # calculates how large each range is that each pollard-rho process can start in
    numRange = math.ceil(math.pow(N, 1/2)//numPollardRho)+1
    lowerBound = 0
    upperBound = numRange

    # starts brute force search
    p1 = Process(target=bruteForce, args=(N, ))
    p1.start()

    # starts first pollard-rho process
    p2 = Process(target=pollardRho, args=(N, lowerBound, upperBound))
    p2.start()

    lowerBound += numRange
    upperBound += numRange

    # starts second pollard-rho process
    p3 = Process(target=pollardRho, args=(N, lowerBound, upperBound))
    p3.start()

    lowerBound += numRange
    upperBound += numRange

    # starts third pollard-rho process
    p4 = Process(target=pollardRho, args=(N, lowerBound, upperBound))
    p4.start()

    lowerBound += numRange
    upperBound += numRange

    # starts fourth pollard-rho process
    p5 = Process(target=pollardRho, args=(N, lowerBound, upperBound))
    p5.start()

    lowerBound += numRange
    upperBound += numRange

    # starts fifth pollard-rho process
    p6 = Process(target=pollardRho, args=(N, lowerBound, upperBound))
    p6.start()

    lowerBound += numRange
    upperBound += numRange

    # starts sixth pollard-rho process
    p7 = Process(target=pollardRho, args=(N, lowerBound, upperBound))
    p7.start()

    lowerBound += numRange
    upperBound += numRange

    # starts seventh pollard-rho process
    p8 = Process(target=pollardRho, args=(N, lowerBound, upperBound))
    p8.start()

    p1.join()
    p2.join()
    p3.join()
    p4.join()
    p5.join()
    p6.join()
    p7.join()
    p8.join()


if __name__ == '__main__':
    print("Enter your desired number to factor: ")
    N = int(input())
    main()
