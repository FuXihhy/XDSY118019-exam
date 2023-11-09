"""
This is the unified lib for Root Finding

Method provided: Bisection, Newton
Applicable function type : scalar
"""

from sys import float_info

class Root:
    """
    This is a simple class to store the calls, status of convergence, error, the numerical root and the possible alarm during solving.
    """
    def __init__(self, method : str = None, converged : bool = None, iterations : int = None, 
                 root : float = None, error : float = None, alarm : list = []) -> None:
        self.method = method
        self.iterations = iterations
        self.converged = converged
        self.root = root
        self.error = error
        self.alarm = alarm
    
    def __str__(self) -> str:
        # format the alarm part
        alarm = ""
        l = len(self.alarm)
        if l == 0:
            alarm = 'None'
        else:
            alarm = self.alarm[0]
            for i in range(1, l):
                alarm += "\n              " + self.alarm[i]
        
        return f"\n\
     method : {self.method}\n\
  converged : {self.converged}\n\
iterracions : {self.iterations}\n\
       root : {self.root}\n\
      error : {self.error}\n\
      alarm : {alarm}\n"

def __raise_alarm(alarm : list = [], Error : str = '') -> str:
    temp_alarm = ""
    for str in alarm:
        temp_alarm += str + '\n'
    raise ValueError(temp_alarm + Error)

def __Bisection(f, a = 0, b = 0, epsilon = 100 * float_info.epsilon,
                 alarm : list = []) -> Root:
    """
    Finding the root using the method of Bisection

    Parameters
    ----------
    f : callable
        must be a scalar funtion to find root of, or there will raise error when calling

    a : float
        the left endpoint of the interval to find root on

    b : float
        the right endpoint of the interval to find root on

    epsilon : float, optional
        the shreshold value or the tolerance of the root, default value is 100 times the sys.float_info.epsilon

    alarm : str, optional
        the alarm passed by the father function to show in the final result
    
    Returns
    -------
    Return the informations inlcuding the convergence, calls of function, error, and the numertical root
    Notice that the approximate root is the average of the last pair of (a, b),
    which indicate that the error = (b - a) / 2 < epsilon
    """

    __sgn = lambda x : 1 if x > 0 else (-1 if x < 0 else 0)
    # use the function `sgn` to skip the operation of mutipling two complex number

    if a > b :
        __raise_alarm(alarm, "Insensible interval")
    a_sgn = __sgn(f(a))
    b_sgn = __sgn(f(b))
    iterations = 0
    mid = 0
    temp = 0
    # difine them out of the loop structure, saving time to create them
    if a_sgn == b_sgn and a_sgn != 0:
        mid = (a + b) / 2
        temp = __sgn(f(mid))
        iterations = 1
        if temp == a_sgn:
            __raise_alarm(alarm, "The values of the function at the endpoints of the bracket aren't of opposite signs! Attempting Failed!")
        else:
            # Attempting to correct the interval with the same sign at the endpoint. Only once!
            b = mid
            b_sgn = __sgn(f(b))
            alarm.append("Given endpoints with same sign! Correcting attempt successed!")

    if a_sgn == 0:
        return Root('Bisection', True, iterations, a, 0, alarm = alarm)
    elif b_sgn == 0:
        return Root('Bisection', True, iterations, b, 0, alarm = alarm)
    
    while (b - a > 2 * epsilon):
        mid = (a + b) / 2
        temp = __sgn(f(mid))
        iterations += 1
        if temp == 0:
            a = mid
            b = mid
            break
        elif temp * a_sgn < 0:
            b = mid
        else:
            a = mid
    return Root('Bisection', True, iterations, mid, (b - a) / 2, alarm = alarm)

def __Newton(f, fprime, x0 = None, epsilon = 100 * float_info.epsilon,
              loop_tol = 1000, alarm : list = []) -> Root:
    """
    Finding the root using the method of Bisection

    Parameters
    ----------
    f : callable
        must be a scalar funtion to find root of, or there will raise error when calling

    fprime : callable
        the derivative of f

    x0 : float
        the initial value to start from

    epsilon : float, optional
        the shreshold value or the tolerance of the root, default value is 100 times the sys.float_info.epsilon

    loop_tol : int, optional
        the tolerance of the numbers of loop when iterating, default as 1000 times
    
    alarm : str, optional
        the alarm passed by the father function to show in the final result
    
    Returns
    -------
    Return the informations inlcuding the convergence, calls of function, error, and the numertical root
    Notice that the approximate root is the average of the last pair of (a, b),
    which indicate that the error = (b - a) / 2 < epsilon
    """
    
    abs = lambda x : x if x > 0 else -x

    if x0 == None:
        alarm.append("No initial value given! Set as 0!")
        x0 = 0
    if not fprime:
        raise ValueError("No derivative given when solving with Newton method!")
    
    iterate_calls = 1
    x1 = None
    temp = f(x0)
    zero_division_flag = False
    if temp == 0:
        x1 = x0
    else:
        temp_prime = fprime(x0)
        if temp_prime == 0:
            # adding a epsilon to prevent the ZeroDivisionError
            if zero_division_flag == False:
                alarm.append("Derivative of f equals zero when iterating! Add epsilon!")
                zero_division_flag = True
            temp_prime += float_info.epsilon
        x1 = x0 - temp / temp_prime
    dic = [x0, x1]
    # We choose to store the points that we've been to, 
    # for the O(n^2) time complexity is sometimes less important compared to some function and the derivative of f that are much more complex to calculate
    while(abs(x0 - x1) > epsilon):
        x0 = x1
        temp = f(x0)
        iterate_calls += 1
        if temp == 0:
            x1 = x0
            return Root('Newton', True, iterations = iterate_calls, root = x1, error = 0, alarm = alarm)
        else:
            temp_prime = fprime(x0)
            if temp_prime == 0:
                # adding a epsilon to prevent the ZeroDivisionError
                if zero_division_flag == False:
                    alarm.append("Derivative of f equals zero when iterating! Add epsilon!")
                    zero_division_flag = True
                temp_prime += float_info.epsilon
            
            x1 = x0 - temp / temp_prime
            # To check if there exist a cycle during solving
            for i in range(iterate_calls - 1):# dic[iterate_calls - 1] = x0
                if abs(dic[i] - x1) <= float_info.epsilon:
                    alarm.append("Cycle exists during solving!")
                    return Root('Newton', False, iterations = iterate_calls, root = 'None', error = 'None', alarm = alarm)
        dic.append(x1)
        if iterate_calls > loop_tol :
            alarm.append("Didn't converge under the limit of numbers of iteration!\nThe root is the latest approximate root while iterating!")
            return Root('Newton', False, iterations = iterate_calls, root = x1, error = 'None', alarm = alarm)
        
    return Root('Newton', True, iterations = iterate_calls, root = dic[iterate_calls], error = (dic[iterate_calls] + dic[iterate_calls - 1]) / 2, alarm=alarm)

def root_finding(f, method :str = None, braket : list = None,
                  fprime = None, x0 = None, x1 = None, 
                  epsilon = 100 * float_info.epsilon, loop_tol = 1000,
                  plot = False, plot_interval = None) -> Root:
    """
    Finding a root of a given scalar function

    Parameters
    ----------
    f : callable
        must be a scalar funtion to find root of, or there will raise ValueError

    method : str, optional
        there is two given method to choose from:
            - Bisection : `bracket` needed, which implies the signum of the endpoints funtion value should be opposite
            - Newton :  `fprime`,`x0` needed; `x1` optional.

    bracket : list of two element, optional
        standing for the begining and the end of a interval(endpoints included)

    fprime : callable, optional
        the derivative of f

    x0 : float, optional
        the default start value to try in the method of `Newton`

    x1 : float, optional
        the second start value to try in the method of `Newton`

    epsilon : float, optional
        the shreshold value or the tolerance of the root, default value is 100 times the sys.float_info.epsilon

    loop_tol : int, optional
        the tolerance of the numbers of the loops while using the `Newton` method

    plot : bool, optinal
        whether to plot the objective function

    plot_interval : list, optional
        the interval to plot on, default as bracket, if given, or [-5, 5]
    """

    if plot:
        if plot_interval == None:
            if braket:
                plot_interval = braket
            else:
                plot_interval = [-5, 5]
            print(f"Ploting interval no given, set as {plot_interval}!\n")
        if plot_interval[0] > plot_interval[1]:
            __raise_alarm([], 'Illegal ploting interval!')
        import matplotlib.pyplot as plt
        import numpy as np
        x = np.arange(plot_interval[0], plot_interval[1] + 0.01, 0.01)
        plt.plot(x, f(x))

    alarm = []
    if method == None:
        if braket:
            method = 'Bisection'
            alarm.append("No method given, automatically chosed as Bisection!")
        elif fprime:
            method = 'Newton'
            alarm.append("No method given, automatically chosed as Newton!")

    if method == 'Bisection':
        return __Bisection(f, braket[0], braket[1], epsilon = epsilon, alarm = alarm)
    
    elif method == 'Newton':
        root = __Newton(f, fprime, x0, epsilon = epsilon, alarm = alarm, loop_tol = loop_tol)
        if x1 and root.converged == False:
            alarm.append(f"First attemption of x_0 = {x0} failed! Here comes the second try with the initial value as {x1}")
            root = __Newton(f, fprime, x1, epsilon = epsilon, alarm = alarm, loop_tol = loop_tol)
        return root
    
    else:
        return Root(alarm=["Illegal operation!"])

if __name__ == "__main__":
    import numpy as np
    def f(x):
        return x ** 2 - 1

    def fprime(x):
        return 2 * x
    #root = root_finding(f, "Bisection", braket=[-1, 1])
    root = root_finding(f, "Newton", fprime=fprime, x0 = 1, plot = True, plot_interval=[-2, 2])
    print(root)