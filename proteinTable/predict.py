import requests
import pandas as pd
import os.path
import xlrd
import sklearn as sk
from sklearn.linear_model import LogisticRegression

data = pd.read_excel('protein_table.xlsx')

data.head()

y = data.iloc[:,34]
print (y)
X = data.iloc[:, [6,8,15,28]] #,30,32]]
print(X)
LR = LogisticRegression(random_state=0, solver='lbfgs', multi_class='ovr').fit(X, y)
print(len(LR.predict(X.iloc[0:,:])))
print(LR.predict(X.iloc[0:,:]))
print(round(LR.score(X,y), 4))

a = pd.DataFrame(LR.predict(X.iloc[0:,:]))
a.to_excel(excel_writer = "predict.xlsx")