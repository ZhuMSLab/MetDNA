"""
created by Xin Xiong<onlybelter@outlook.com>, https://github.com/OnlyBelter
"""

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os


ERROR_CODE = {0: 'valid',
              1: '[ERROR] - ' + 'Has missing value or whitespace',
              2: '[ERROR] - ' + 'The first column name of MS1 Data is not "name"',
              3: '[ERROR] - ' + 'The second column name of MS1 Data is not "mz"',
              4: '[ERROR] - ' + 'The third column name of MS1 Data is not "rt"',
              5: '[WARNING] - ' + 'max(rt) <= 60, rt\'s unit may be not "second"',
              6: '[WARNING] - ' + 'Some peaks\' have more than 50% of zero values in all samples',
              7: '[ERROR] - ' + 'The first column name of Sample Infomation is not "sample.name"',
              8: '[ERROR] - ' + 'The second column name of Sample Infomation is not "group"',
              9: '[ERROR] - ' + 'Sample names in MS Data and Sample Information are different'}


def check_data(file_type, file_handle, save_path, pol):
    """
    only check ms1_data and sample_info
    :param file_type: ms1_data, sample_info
    :param file_handle: file comes from user web input
    :param save_path: save 'peak intensity profile' file
    :param pol: polarity, neg or pos
    :return:
    """
    error_code = []
    file_info = {}
    valid = True
    # if file_type == 'ms1_data' or file_type == 'sample_info':
    if file_type == 'ms1_data' or file_type == 'sample_info':
        file_data = pd.read_csv(file_handle)
        column_name = file_data.columns.tolist()
        file_data = file_data.replace(r'\s+', np.nan, regex=True)
        has_nan = file_data.isnull().values.any()
        if has_nan:
            error_code.append(1)
            valid = False
        if file_type == 'ms1_data':
            if column_name[0] != 'name':
                error_code.append(2)
                valid = False
            if column_name[1] != 'mz':
                error_code.append(3)
                valid = False
            if column_name[2] != 'rt':
                error_code.append(4)
                valid = False

            if valid:
                if file_data.rt.max() <= 60:
                    error_code.append(5)
                sample_int = file_data.iloc[:, 3:]
                sample_size = len(sample_int.columns.tolist())
                zero_percent = ((sample_int == 0).sum(axis=1))/float(sample_size)
                col_num = (zero_percent >= 0.5).sum()  # column number that the column has more than 50% 0 value
                if col_num > 0:
                    error_code.append(6)
                file_info['sample_size'] = sample_size
                file_info['sample_name'] = column_name[3:]
                file_info['peak_size'] = file_data.shape[0]
                if not os.path.exists(save_path):
                    os.makedirs(save_path)
                path = os.path.join(save_path, 'peak_intensity_profile_{}.png'.format(pol))
                plot_peaks(file_data, path, pol)
                # file_info['peak_profile_path'] = path

        if file_type == 'sample_info':
            if column_name[0] != 'sample.name':
                error_code.append(7)
                valid = False
            if column_name[1] != 'group':
                error_code.append(8)
                valid = False
            if valid:
                # no any error, store file info
                groups = dict(list(file_data.groupby('group')))
                for key in groups:
                    _key = key
                    if type(key) is not str:
                        try:
                            _key = str(key)
                        except:
                            _key = ''
                    if _key:
                        file_info[_key] = groups[key]['sample.name'].tolist()

        if valid:
            error_code.append(0)  # data is valid

    return {'error_code': error_code, 'file_info': file_info}


def plot_peaks(df, save_path, pos):
    # https://stackoverflow.com/a/47074245/2803344
    # plt = matplotlib.pyplot
    _ = pd.DataFrame(df.iloc[:, :])  # df is a dataFrame
    _['int_mean'] = _.iloc[:, 3:].mean(axis=1)
    _['int_mean_log2'] = np.log2(_['int_mean'])
    fig, ax = plt.subplots()
    fig.subplots_adjust(hspace=0.5, left=0.12, right=0.93)
    hb = ax.hexbin(x=_['rt'], y=_['mz'], C=_['int_mean_log2'],
                   reduce_C_function=np.max, gridsize=100)
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label('log2(Mean of Peak Intensity)')
    plt.xlabel('Retention Time(RT, second)')
    plt.ylabel('Mass-to-charge Ratio(m/z)')
    plt.title('Peak Intensity Profile({}.)'.format(pos))
    # save_path = os.path.join(save_path, )
    # if not os.path.exists(save_path):
    plt.savefig(save_path, dpi=200)


# dir = r'E:\rd\014_metDNA\demo\MetDNA.demo.data.NEG\NEG'
# # file_name = 'sample.info.csv'
# file_name = 'data.csv'
# file_path = os.path.join(dir, file_name)
# # file_type = 'sample_info'
# file_type = 'ms1_data'
# check_result = check_data(file_type, file_path)
# print(check_result)
# plot_peaks(file_path)
