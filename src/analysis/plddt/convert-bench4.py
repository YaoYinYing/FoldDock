#!/usr/bin/env python3
import numpy as np
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#import seaborn as sns
#import re
import pandas as pd

df=pd.read_csv("summary.csv",sep=",",index_col=False)

df["id1"]=df.name.str[:4]+"_A"
df["id2"]=df.name.str[:4]+"_B"
df1=df.loc[df["rank"]==1][["id1","id2","pcd","mmd"]]
df2=df.loc[df["rank"]==2][["id1","id2","pcd","mmd"]]
df3=df.loc[df["rank"]==3][["id1","id2","pcd","mmd"]]
df4=df.loc[df["rank"]==4][["id1","id2","pcd","mmd"]]
df5=df.loc[df["rank"]==5][["id1","id2","pcd","mmd"]]
df1=df1.merge(df2,on=["id1","id2"],suffixes=("_1","_2"))
df1=df1.merge(df3,on=["id1","id2"],suffixes=("_2","_3"))
df1=df1.merge(df4,on=["id1","id2"],suffixes=("_3","_4"))
df1=df1.merge(df5,on=["id1","id2"],suffixes=("_4","_5"))
df1=df1.rename(columns={"pcd":"pcd_5","mmd":"mmd_5"})

df1.to_csv("pconsdock-bench4.csv")
