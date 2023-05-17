import config
import numpy as np
import matplotlib.pylab as plt


def main():
    t_step = 1 / 100
    t_start = 0
    t_end = 12 * 15

    # N = 3_444_000_000 # TODO: mis wel N nemen op N van 2022
    # # N = 5 * 10 ** 9
    # S0 = 1_144_000_000
    # I0 = 2_000_000_000
    # R0 = 300_000_000
    # alpha = 0.99
    # beta = 0.32
    # # beta = 0.02
    N = 3_444 * 10**6
    I0 = 1_936 * 10**6
    R0 = 0.15 * I0
    S0 = N - I0 - R0
    downloads_per_month = (660 * 10**6) / 12 # \frac{d}{\Delta t}
    average_infection_rate = (2_129 * 10**6 - 1_936 * 10**6) / 12 # \frac{\Delta u}{\Delta t}
    alpha = N / (S0 * I0) * downloads_per_month
    beta = (alpha * S0) / N - average_infection_rate / I0

    n_end = int(1 + (t_end - t_start) / t_step) # FIXME:
    # error if (t_end - t_start) // t_step =/= int
    t_arr = np.linspace(t_start, t_end, n_end, endpoint=True)

    S_arr = np.zeros(n_end)
    I_arr = np.zeros(n_end)
    R_arr = np.zeros(n_end)

    S_arr[t_start] = S0
    I_arr[t_start] = I0
    R_arr[t_start] = R0

    for n, t_n in enumerate(t_arr[:-1]):
        S_arr[n + 1] = S_arr[n] + t_step * (-alpha / N * S_arr[n] * I_arr[n])
        I_arr[n + 1] = I_arr[n] + t_step * (alpha / N * S_arr[n] * I_arr[n] - beta * I_arr[n])
        R_arr[n + 1] = R_arr[n] + t_step * (beta * I_arr[n])

    plt.plot(t_arr, S_arr, color='blue')
    plt.plot(t_arr, I_arr, color='orange')
    plt.plot(t_arr, R_arr, color='green')
    plt.show()


if __name__ == '__main__':
    main()