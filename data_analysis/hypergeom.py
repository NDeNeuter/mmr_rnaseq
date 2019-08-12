#!/usr/bin/env python

from math import log

def factorial(n):
    
    ''' Calculates factorial n = n*(n-1)*(n-2)*...*1. '''
    
    out = float(1)
    for i in range(1,n+1):
        out = out*i
        
    return out


def partial_factorial(n, z):
    
    ''' Calculates a partial factorial for n = n*(n-1)*...*(z+1). If z = 0, this becomes a basic factorial again. '''

    out = float(1)
    for i in range(z+1, n+1):
        out = out*i
        
    return out


def draw(x, N):
    
    ''' Calculates the number of possible ways to draw x items from a population of N items without replacement. '''
    
    n = N
    k = x
    if k >= n-k:
        out = partial_factorial(n, k)/factorial(n-k)
    else:
        out = partial_factorial(n, n-k)/factorial(k)
    
    return out


def pmf(s_sample, N_sample, s_pop, N_pop):

    ''' Probability of drawing exactly s_sample successes from an N_sample sample population given s_pop successes in a total population of size N_pop. '''

    if s_sample > s_pop or N_sample-s_sample > N_pop-s_pop:
        probability = 0
    else:
        succes = draw_log(s_sample, s_pop)
        fail = draw_log(N_sample-s_sample, N_pop-s_pop)
        total = draw_log(N_sample, N_pop)
        probability = float(succes)+fail-total

        return 10**probability

    
def cdf_less(s_sample, N_sample, s_pop, N_pop):

    ''' Probability of drawing s_sample successes or less from an N_sample sample population given s_pop successes in a total population of size N_pop. '''

    s_sample = int(s_sample)
    N_sample = int(N_sample)
    s_pop = int(s_pop)
    N_pop = int(N_pop)
    
    probability = 0
    for i in range(s_sample+1):
        pmf_i = pmf(i, N_sample, s_pop, N_pop)
        if pmf_i != None:
            probability += pmf_i

    return probability

    
def cdf_more(s_sample, N_sample, s_pop, N_pop):

    ''' Probability of drawing s_sample successes or more from an N_sample sample population given s_pop successes in a total population of size N_pop. '''

    s_sample = int(s_sample)
    N_sample = int(N_sample)
    s_pop = int(s_pop)
    N_pop = int(N_pop)
    
    probability = 0
    for i in range(s_sample, N_sample+1):
        pmf_i = pmf(i, N_sample, s_pop, N_pop)
        if pmf_i != None:
            probability += pmf_i

    return probability


def factorial_log(n):
    
    ''' Calculates factorial n = n*(n-1)*(n-2)*...*1.
    To avoid large numbers, returns the log of the factorial. '''
    
    out = float(0)
    for i in range(1,n+1):
        out += log(i, 10)
        
    return out
    
    
def partial_factorial_log(n, z):
    
    ''' Calculates a partial factorial for n = n*(n-1)*...*(z+1). If z = 0, this becomes a basic factorial again.
    To avoid large numbers, returns the log of the partial factorial. '''

    out = float(0)
    for i in range(z+1, n+1):
        out += log(i, 10)
        
    return out


def draw_log(x, N):
    
    ''' Calculates the number of possible ways to draw x items from a population of N items without replacement.
    To avoid large numbers, returns the log of the number of possible draws. '''
    
    n = N
    k = x
    if k >= n-k:
        out = partial_factorial_log(n, k)-factorial_log(n-k)
    else:
        out = partial_factorial_log(n, n-k)-factorial_log(k)
    
    return out
