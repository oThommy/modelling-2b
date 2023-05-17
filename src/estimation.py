def main():
    # TODO: move to config?
    N = 3_444 * 10**6
    I0 = 1_936 * 10**6
    R0 = 0.15 * I0
    S0 = N - I0 - R0
    downloads_per_month = (660 * 10**6) / 12 # \frac{d}{\Delta t}
    average_infection_rate = (2_129 * 10**6 - 1_936 * 10**6) / 12 # \frac{\Delta u}{\Delta t}
    alpha = N / (S0 * I0) * downloads_per_month
    beta = (alpha * S0) / N - average_infection_rate / I0
    
    print(f'N={N / 1e6}+E06')
    print(f'I0={I0 / 1e6}+E06')
    print(f'R0={R0 / 1e6}+E06')
    print(f'S0={S0 / 1e6}+E06')
    print(f'downloads_per_month={downloads_per_month / 1e6}+E06')
    print(f'average_infection_rate={average_infection_rate / 1e6}+E06')
    print(f'{alpha=}')
    print(f'{beta=}')


if __name__ == '__main__':
    main()