from pathlib import Path
from dataclasses import dataclass
from operator import attrgetter
import datetime as dt
import utils


ROOT_DIR = Path(__file__).parent.parent
SRC_DIR = ROOT_DIR / 'src'
OUT_DIR = ROOT_DIR / 'out'
FIGURES_DIR = OUT_DIR / 'figures'
DATA_DIR = ROOT_DIR / 'data'

for dir in [OUT_DIR, FIGURES_DIR, DATA_DIR]:
    utils.ensure_dir_exists(dir)

ENCODING = 'utf-8'

DAYS_IN_MEAN_GREGORIAN_YEAR = 365.2425
MILLISECONDS_IN_MONTH = DAYS_IN_MEAN_GREGORIAN_YEAR * 24 * 60 * 60 * 1000 / 12

def __get_estimates_facebook_sir_constants():
    N = 3_444 * 10**6
    I0 = 1_936 * 10**6
    R0 = 0.15 * I0
    S0 = N - I0 - R0
    downloads_per_month = (660 * 10**6) / 12 # \frac{D}{\Delta t}
    average_infection_rate = (2_129 * 10**6 - 1_936 * 10**6) / 12 # \frac{\Delta I}{\Delta t}
    alpha = N / (S0 * I0) * downloads_per_month
    beta = (alpha * S0) / N - average_infection_rate / I0
    t_start_milliseconds = dt.datetime(2017, 1, 1).timestamp() * 1000

    return N, I0, R0, S0, downloads_per_month, average_infection_rate, alpha, beta, t_start_milliseconds

@dataclass(frozen=True)
class FacebookSirConstants:
    N: float
    I0: float
    R0: float
    S0: float
    downloads_per_month: float
    average_infection_rate: float
    alpha: float
    beta: float
    t_start_milliseconds: float

FACEBOOK_SIR_CONSTANTS = FacebookSirConstants(*__get_estimates_facebook_sir_constants())

def main():
    N, I0, R0, S0, downloads_per_month, average_infection_rate, alpha, beta = attrgetter('N', 'I0', 'R0', 'S0', 'downloads_per_month', 'average_infection_rate', 'alpha', 'beta')(FACEBOOK_SIR_CONSTANTS)

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