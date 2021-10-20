import pandas as pd
import numpy as np
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LinearRegression
from sklearn.metrics import (
    r2_score,
    mean_absolute_error,
)
from scipy.stats import pearsonr


degree=2

# dataset = pd.read_csv('prositIRTs.csv')
# X = np.array(dataset['prositIRT']).reshape(-1, 1)
# y = dataset['rt']



# polyreg.fit(X,y)



KEYS = {
    'MeV': [
        'MeV_BLCL_allFractions.sample{}.csv',
        'MeV_GRLCL_allFractions.sample{}.csv',
        'MeV_MA0009-BE08_allFractions.sample{}.csv',
        'MeV_MA0009-BE09_allFractions.sample{}.csv',
        'MeV_MB211_KCL_allFractions.sample{}.csv',
    ],
    'SPL_90-10': [
        '20200228.SPL1.sample{}.csv',
        '20200228.SPL2.sample{}.csv',
        '2020316.SPL3.sample{}.csv',
        '2020316.SPL4.sample{}.csv',
        '2020316.SPL5.sample{}.csv',
        '2020316.SPL6.sample{}.csv',
        '2020316.SPL7.sample{}.csv',
        '2020316.SPL8.sample{}.csv',
    ],
    'SPL_500': [
        '20200228.SPL1.sample{}.csv',
        '20200228.SPL2.sample{}.csv',
        '2020316.SPL3.sample{}.csv',
        '2020316.SPL4.sample{}.csv',
        '2020316.SPL5.sample{}.csv',
        '2020316.SPL6.sample{}.csv',
        '2020316.SPL7.sample{}.csv',
        '2020316.SPL8.sample{}.csv',
    ],
}
SEQUENCE_KEY = 'Mascot_seq'
RT_KEY = 'RT'

pearson_vals = []
r2_vals = []
mae_vals = []
df_keys = []
source_keys = []

base_folder = 'gpuProsit_RT_benchmark'
for group in ['MeV', 'SPL_90-10', 'SPL_500']:
    sub_folder = base_folder + '/' + group
    sub_folder_2 = sub_folder + '/{}'

    for idx in range(1, 11):
        for div in KEYS[group]:

            train_data_key = sub_folder_2.format('train') + '/prositPred' + div.format(idx)
            test_data_key = sub_folder_2.format('test') + '/prositPred' + div.format(idx)

            train_dataset = pd.read_csv(train_data_key)
            test_dataset = pd.read_csv(test_data_key)

            # drop_lamdba = lambda x : False if 'C' in x else True
            # train_drop = train_dataset['Mascot_seq'].apply(drop_lamdba)
            # train_dataset = train_dataset[train_drop]
            # test_drop = test_dataset['Mascot_seq'].apply(drop_lamdba)
            # test_dataset = test_dataset[test_drop]

            x_train = np.array(train_dataset['prositIRT']).reshape(-1, 1)
            y_train = train_dataset[RT_KEY]

            # polyreg=LinearRegression()
            polyreg=make_pipeline(PolynomialFeatures(degree),LinearRegression())
            polyreg.fit(x_train, y_train)

            x_test = np.array(test_dataset['prositIRT']).reshape(-1, 1)
            y_test = test_dataset[RT_KEY]

            pred_y = polyreg.predict(x_test)

            pearson_value = pearsonr(y_test, pred_y)
            r2_value = r2_score(y_test, pred_y)
            mae = mean_absolute_error(y_test, pred_y)
            pearson_vals.append(pearson_value[0])
            r2_vals.append(r2_value)
            mae_vals.append(mae)
            source_keys.append(group)
            df_keys.append(div.format(idx))




            # import matplotlib.pyplot as plt

            # plt.figure()
            # plt.scatter(x_test, y_test)
            # plt.plot(x_test, polyreg.predict(x_test),color="black")
            # plt.title("Polynomial regression with degree "+str(degree))
            # plt.show()
# print(len(source_keys))
all_results = pd.DataFrame(
    {
        'source': pd.Series(source_keys),
        'df_key': pd.Series(df_keys),
        'pearson': pd.Series(pearson_vals),
        'r2': pd.Series(r2_vals),
        'mae': pd.Series(mae_vals)
    }        
)
print(all_results.median())
print(all_results.mean())
# print(all_results.max())
# print(all_results.min())
print(all_results.groupby('source').mean())
all_results.to_csv('gpuProsit_RT_benchmark/summaryResults.csv', index=False)