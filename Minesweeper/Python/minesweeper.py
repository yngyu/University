import numpy as np
import copy
import random
import queue

def first_input_linenumber():
    while True:
        print("縦のマスの数を入れてください。大きすぎると上手く表示出来なくなります。")
        temp = input("")
        if temp.isdecimal() and temp[0] != "0":
            return int(temp)
        else:
            print("自然数を入れてください")
            continue

def first_input_columnnumber():
    while True:
        print("横のマスの数を入れてください。大きすぎると上手く表示出来なくなります。")
        temp = input("")
        if temp.isdecimal() and temp[0] != "0":
            return int(temp)
        else:
            print("自然数を入れてください")
            continue
def first_input_bombnumber(linenumber,columnnumber):
    while True:
        print("地雷の数を入れてください")
        temp = input("")
        if temp.isdecimal() and temp[0] != "0":
            if int(temp) >= linenumber*columnnumber-8:
                print("地雷の数が多すぎます")
                continue
            else:
                return int(temp)
        else:
            print("自然数を入れてください")
            continue

linenumber = first_input_linenumber()
columnnumber = first_input_columnnumber()
bombnumber = first_input_bombnumber(linenumber,columnnumber)
display_board_array = np.full((linenumber,columnnumber),"o")
data_board_array = np.zeros((linenumber,columnnumber))

def know_digit(number):
    digit = 0
    while number:
        number //= 10
        digit += 1
    return digit

def display_board_func(display_board_array,linenumber,columnnumber):
    line_digit = know_digit(linenumber)
    column_digit = know_digit(columnnumber)
    for line in range(linenumber+1):
        if line == 0:
            for i in range(line_digit):
                print(" ",end="")
            for column in range(1,columnnumber+1):
                print(column,end="")
                for i in range(column_digit-know_digit(column)):
                    print(" ",end="")
            print()
        else:
            print(line,end="")
            for i in range(line_digit-know_digit(line)):
                print(" ",end="")
            for column in range(1,columnnumber+1):
                print(display_board_array[line-1][column-1],end="")
                for i in range(column_digit-1):
                    print(" ",end="")
            print()

display_board_func(display_board_array,linenumber,columnnumber)

def input_linenumber(linenumber):
    while True:
        print("マスの位置の縦の座標を入れてください")
        temp = input("")
        if temp.isdecimal() and temp[0] != "0":
            if int(temp) <= linenumber:
                return int(temp)
            else:
                print("そのマスはありません")
                continue
        else:
            print("自然数を入れてください")
            continue

def input_columnnumber(columnnumber):
    while True:
        print("マスの位置の横の座標を入れてください")
        temp = input("")
        if temp.isdecimal() and temp[0] != "0":
            if int(temp) <= columnnumber:
                return int(temp)
            else:
                print("そのマスはありません")
                continue
        else:
            print("自然数を入れてください")
            continue

def count_around_bomb(data_board_array,linenumber,columnnumber):
    threenumber = [-1,0,1]
    for line in range(linenumber):
        for column in range(columnnumber):
            if data_board_array[line][column] == -1:
                continue
            else:
                for i in threenumber:
                    for j in threenumber:
                        if i == 0 and j == 0:
                            continue
                        elif line+i < 0 or column+j < 0:
                            continue
                        else:
                            try:
                                if data_board_array[line+i][column+j] == -1:
                                    data_board_array[line][column] += 1
                            except IndexError:
                                continue

def find_blank_space(display_board_array,data_board_array,line_dig_point,column_dig_point):
    find_blank = queue.Queue()
    find_blank.put([line_dig_point,column_dig_point])
    while not find_blank.empty():
        temp = find_blank.get()
        line = temp[0]
        column = temp[1]
        threenumber = [-1,0,1]
        for i in threenumber:
            for j in threenumber:
                if line+i < 0 or column+j < 0:
                    continue
                else:
                    try:
                        if display_board_array[line+i][column+j] == "o":
                            if data_board_array[line+i][column+j] == 0:
                                display_board_array[line+i][column+j] = " "
                                find_blank.put([line+i,column+j])
                            else:
                                display_board_array[line+i][column+j] = str(data_board_array[line+i][column+j])
                        else:
                            continue
                    except IndexError:
                        continue


def initialize_board_func(display_board_array,data_board_array,linenumber,columnnumber,bombnumber):
    print("最初に掘るマスの場所を入れてください")
    line_dig_point = input_linenumber(linenumber)
    line_dig_point -= 1
    column_dig_point = input_columnnumber(columnnumber)
    column_dig_point -= 1
    display_board_array[line_dig_point][column_dig_point] = " "
    decide_where_is_bomb = []
    around_the_first_dig_point = []
    threenumber = [-1,0,1]
    for i in threenumber:
        for j in threenumber:
            if line_dig_point+i < 0 or column_dig_point+j < 0:
                continue
            try:
                if data_board_array[line_dig_point+i][column_dig_point+j] == 0:
                    around_the_first_dig_point.append((line_dig_point+i)*columnnumber+column_dig_point+j)
            except IndexError:
                continue
    for i in range(linenumber*columnnumber):
        if i in around_the_first_dig_point:
            continue
        else:
            decide_where_is_bomb.append(i)
    random.shuffle(decide_where_is_bomb)
    for i in range(bombnumber):
        data_board_array[decide_where_is_bomb[i]//columnnumber][decide_where_is_bomb[i]%columnnumber] = -1
    count_around_bomb(data_board_array,linenumber,columnnumber)
    find_blank_space(display_board_array,data_board_array,line_dig_point,column_dig_point)

initialize_board_func(display_board_array,data_board_array,linenumber,columnnumber,bombnumber)

def judge_game_clear(display_board_array,data_board_array,linenumber,columnnumber):
    for line in range(linenumber):
        for column in range(columnnumber):
            if data_board_array[line][column] != -1:
                if display_board_array[line][column] == "o" or display_board_array[line][column] == "F":
                    return False
                else:
                    continue
            else:
                continue
    return True

def display_bomb(display_board_array,data_board_array,linenumber,columnnumber):
    for line in range(linenumber):
        for column in range(columnnumber):
            if data_board_array[line][column] == -1:
                display_board_array[line][column] = "x"
            else:
                continue

sparebombnumber = copy.deepcopy(bombnumber)
while True:
    if judge_game_clear(display_board_array,data_board_array,linenumber,columnnumber):
        display_board_func(display_board_array,linenumber,columnnumber)
        print("おめでとうございます！ゲームクリアです！")
        break
    display_board_func(display_board_array,linenumber,columnnumber)
    print("残りの地雷数:",end="")
    print(sparebombnumber)
    while True:
        print("d:掘る, f;旗を立てる, c:旗を消す")
        command = input("行動を入力してください:")
        if command == "d":
            print("どこを掘りますか?")
            break
        elif command == "f":
            print("どこに旗を建てますか?")
            break
        elif command == "c":
            print("どこの旗を消しますか?")
            break
        else:
            print("入力が間違っています")
            continue
    line_dig_point = input_linenumber(linenumber)
    line_dig_point -= 1
    column_dig_point = input_columnnumber(columnnumber)
    column_dig_point -= 1
    if command == "d":
        if display_board_array[line_dig_point][column_dig_point] == "F":
            print("旗があります！")
            continue
        elif display_board_array[line_dig_point][column_dig_point] != "o":
            print("そこは掘れません")
            continue
        elif data_board_array[line_dig_point][column_dig_point] == -1:
            display_board_array[line_dig_point][column_dig_point] = "x"
            display_board_func(display_board_array,linenumber,columnnumber)
            print("残念！地雷がありました...")
            display_bomb(display_board_array,data_board_array,linenumber,columnnumber)
            display_board_func(display_board_array,linenumber,columnnumber)
            print("GAME OVER")
            break
        elif data_board_array[line_dig_point][column_dig_point] == 0:
            find_blank_space(display_board_array,data_board_array,line_dig_point,column_dig_point)
        else:
            display_board_array[line_dig_point][column_dig_point] = str(data_board_array[line_dig_point][column_dig_point])
    elif command == "f":
        if display_board_array[line_dig_point][column_dig_point] != "o":
            print("そこに旗は立てられません")
            continue
        else:
            display_board_array[line_dig_point][column_dig_point] = "F"
            sparebombnumber -= 1
    elif command == "c":
        if display_board_array[line_dig_point][column_dig_point] != "F":
            print("そこに旗はありません")
            continue
        else:
            display_board_array[line_dig_point][column_dig_point] = "o"
            sparebombnumber += 1
    else:
        print("入力が間違っています")