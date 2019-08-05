import numpy as np

def input_linenumber():
    while True:
        print("横のマスの数を入れてください。大きすぎると上手く表示出来なくなります。")
        temp = input("")
        if temp.isdecimal() and temp[0] != "0":
            return temp
        else:
            print("自然数を入れてください")
            continue

def input_columnnumber():
    while True:
        print("縦のマスの数を入れてください。大きすぎると上手く表示出来なくなります。")
        temp = input("")
        if temp.isdecimal() and temp[0] != "0":
            return temp
        else:
            print("自然数を入れてください")
            continue
def input_bombnumber(linenumber,columnnumber):
    while True:
        print("縦のマスの数を入れてください。大きすぎると上手く表示出来なくなります。")
        temp = input("")
        if temp.isdecimal() and temp[0] != "0":
            if temp >= 
        else:
            print("自然数を入れてください")
            continue
linenumber = input_linenumber()
columnnumber = input_columnnumber()
bombnumber = input_bombnumber(linenumber,columnnumber)