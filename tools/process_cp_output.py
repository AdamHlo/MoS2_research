import pandas as pd
import matplotlib.pyplot as plt


def process_evp(evp_path) -> pd.DataFrame:
    with open(evp_path) as file:
        lines = file.readlines()

    columns = lines[0].split()[1:]

    data = []
    for line in lines[1:]:
        split = line.split()
        row = [int(split[0]) if i == 0 else float(split[i]) if i < len(split) else float('NaN') for i in
               range(len(columns))]
        data.append(row)

    dataframe = pd.DataFrame.from_records(data, columns=columns)

    return dataframe


def process_lattice_seq(cel_path, pos_path, vel_path):
    pass


def plot_run():
    fig1 = create_fig('res/4x4_s_vacancy_randomized_500K/mos2.evp')
    fig2 = create_fig('res/4x4_s_vacancy_randomized_500K_additional_scaling/mos2.evp')
    fig3 = create_fig('res/4x4_s_vacancy_randomized_500K_longer_scaling/mos2.evp')
    plt.show()


def create_fig(evp_path):
    evp = process_evp(evp_path)
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    ax1.plot('time(ps)', 'ekinc', data=evp)
    ax1.set_xlabel('time [ps]')
    ax1.set_ylabel('Ry')
    ax1.legend()
    ax2.plot('time(ps)', 'Tion(K)', data=evp)
    ax2.set_xlabel('time [ps]')
    ax2.set_ylabel('K')
    ax2.legend()
    ax3.plot('time(ps)', 'etot', data=evp)
    ax3.plot('time(ps)', 'econt', data=evp)
    ax3.plot('time(ps)', 'econs', data=evp)
    ax3.set_xlabel('time [ps]')
    ax3.set_ylabel('Ry')
    ax3.legend()
    ax4.plot('time(ps)', 'enthal', data=evp)
    ax4.set_xlabel('time [ps]')
    ax4.set_ylabel('Ry')
    ax4.legend()
    fig.tight_layout()
    return fig


def main():
    plot_run()


if __name__ == '__main__':
    main()
