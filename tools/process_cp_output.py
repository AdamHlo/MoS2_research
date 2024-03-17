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


def main():
    evp_path = 'res/mos2.evp'
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
    plt.show()


if __name__ == '__main__':
    main()
