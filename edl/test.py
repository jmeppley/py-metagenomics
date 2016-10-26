import logging
import sys


def myAssertEq(a, b):
    myAssert((a == b), "%s is not equal to %s" % (str(a), str(b)))


def myAssertIs(a, b):
    myAssert((a is b), "%s is not %s" % (str(a), str(b)))


def myAssert(test, msg):
    if not test:
        sys.stderr.write(msg + "\n")
        raise AssertionError
