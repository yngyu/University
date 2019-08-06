#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <queue>
#include <utility>

using namespace std;

int input_linenumber();
int input_columnnumber();
int input_bombnumber();

int know_digit(int); //桁を知る。

void display_board(vector<vector<char>> &board, int linenumber, int columnnumber);

void initialize_board(vector<vector<char>> &display_board_vector, vector<vector<int>> &data_board_vector, int linenumber, int columnnumber, int bombnumber);
void count_around_bomb(vector<vector<int>> &data_board_vector, int linenumber, int columnnumber);

void find_blank_space(vector<vector<char>> &display_board_vector, vector<vector<int>> &data_board_vector, int line_dig_point, int column_dig_point);

bool judge_game_clear(vector<vector<char>> &display_board_vector, vector<vector<int>> &data_board_vector, int linenumber, int columnnumber);

int main()
{
    const int linenumber = input_linenumber();
    const int columnnumber = input_columnnumber();
    const int bombnumber = input_bombnumber();

    int spare_bomb_number = bombnumber;
    
    vector<vector<char>> display_board_vector(linenumber,vector<char>(columnnumber,'o'));
    //旗はF、空は あるのは数字、まだ掘ってないのはoにする。
    vector<vector<int>> data_board_vector(linenumber,vector<int>(columnnumber,0));

    display_board(display_board_vector, linenumber, columnnumber);
    initialize_board(display_board_vector, data_board_vector, linenumber, columnnumber, bombnumber);

    char command;
    int line_dig_point,column_dig_point;
    while(true){
        if(judge_game_clear(display_board_vector, data_board_vector, linenumber,columnnumber)){
            display_board(display_board_vector, linenumber, columnnumber);
            cout << "Conguratulation! Clear!" << endl;
            break;
        }
        display_board(display_board_vector, linenumber, columnnumber);
        cout << "sparebombnumber:";
        cout << spare_bomb_number << endl;
        while(true){
            cout << "d:dig, f:bulid a flag, c:cancel a flag" << endl;
            cout << "input command" << endl;
            cin >> command;
            if(command == 'd'){
                cout << "where to dig" << endl;
                break;
            }else if(command == 'f'){
                cout << "where to build a flag" << endl;
                break;
            }else if(command == 'c'){
                cout << "which flag cancel" << endl;
            }else{
                cout << "input is incorrect" << endl;
            }
        }
        while(true){
            line_dig_point = input_linenumber();
            --line_dig_point;
            column_dig_point = input_columnnumber();
            --column_dig_point;
            try{
                if(display_board_vector.at(line_dig_point).at(column_dig_point)){
                    break;
                }
            }catch(exception &e){
                cout << "input is incorrect" << endl;
                continue;
            }
        }
        if(command == 'd'){
            if(display_board_vector.at(line_dig_point).at(column_dig_point) == 'F'){
                cout << "there is a flag!" << endl;
                continue;
            }else if(display_board_vector.at(line_dig_point).at(column_dig_point) != 'o'){
                cout << "there can't dig" << endl;
                continue;
            }else if(data_board_vector.at(line_dig_point).at(column_dig_point) == -1){
                display_board_vector.at(line_dig_point).at(column_dig_point) = 'x';
                display_board(display_board_vector, linenumber, columnnumber);
                cout << "there is bomb..." << endl;
                cout << "Game Over" << endl;
                break;
            }else if(data_board_vector.at(line_dig_point).at(column_dig_point) == 0){
                find_blank_space(display_board_vector, data_board_vector, line_dig_point, column_dig_point);
            }else{
                display_board_vector.at(line_dig_point).at(column_dig_point) = '0' + data_board_vector.at(line_dig_point).at(column_dig_point);
            }
        }else if(command == 'f'){
            if(display_board_vector.at(line_dig_point).at(column_dig_point) != 'o'){
                cout << "there can't build a flag" << endl;
                continue;
            }else{
                display_board_vector.at(line_dig_point).at(column_dig_point) = 'F';
                --spare_bomb_number;
            }
        }else if(command == 'c'){
            if(display_board_vector.at(line_dig_point).at(column_dig_point) != 'F'){
                cout << "there isn't a flag" << endl;
                continue;
            }else{
                display_board_vector.at(line_dig_point).at(column_dig_point) = 'o';
                ++spare_bomb_number;
            }
        }else{
            cout << "input is incorrect" << endl;
            continue;
        }
    }
    return 0;
}

int input_linenumber()
{
    int linenumber;
    cout << "input number of line" << endl;
    cin >> linenumber;
    return linenumber;
}

int input_columnnumber()
{
    int columnnumber;
    cout << "input number of column" << endl;
    cin >> columnnumber;
    return columnnumber;
}

int input_bombnumber()
{
    int bombnumber;
    cout << "input number of bomn" << endl;
    cin >> bombnumber;
    return bombnumber;
}

void display_board(vector<vector<char>> &board, int linenumber, int columnnumber)
{
    const int line_digit = know_digit(linenumber);
    const int column_digit = know_digit(columnnumber);
    int line,column,i;
    for(line = 0;line <= linenumber;++line){
        if(line == 0){
            for(i = 0;i < line_digit;++i){
                cout << ' ';
            }
            for(column = 1;column <= columnnumber;++column){
                cout << column;
                for(i = 0;i < column_digit - know_digit(column);++i){
                    cout << ' ';
                }
            }
            cout << endl;
        }else{
            cout << line;
            for(i = 0;i < line_digit - know_digit(line);++i){
                cout << ' ';
            }
            for(column = 1;column <= columnnumber;++column){
                cout << board.at(line-1).at(column-1);
                for(i = 0;i < column_digit - 1;++i){
                    cout << ' ';
                }
            }
            cout << endl;
        }
    }
}

void initialize_board(vector<vector<char>> &display_board_vector, vector<vector<int>> &data_board_vector, int linenumber, int columnnumber, int bombnumber)
{
    int line_dig_point,column_dig_point;
    cout << "input first dig point" << endl;
    while(true){
        line_dig_point = input_linenumber();
        --line_dig_point;
        column_dig_point = input_columnnumber();
        --column_dig_point;
        try{
            if(display_board_vector.at(line_dig_point).at(column_dig_point)){
                break;
            }
        }catch(exception &e){
            cout << "input is incorrect" << endl;
            continue;
        }
    }
    display_board_vector.at(line_dig_point).at(column_dig_point) = ' ';
    //爆弾を-1で表すことにする。
    vector<int> decide_where_is_bomb;
    int i,j;
    vector<int> around_the_first_dig_point;
    for(i = -1;i <=1;++i){
        for(j = -1;j <= 1;++j){
            try{
                if(data_board_vector.at(line_dig_point+i).at(column_dig_point+j) == 0){
                    around_the_first_dig_point.push_back((line_dig_point+i)*columnnumber+column_dig_point+j);
                }
            }catch(exception& e){
                continue;
            }
        }
    }
    for(i = 0;i <= linenumber*columnnumber-1;++i){
        if(binary_search(around_the_first_dig_point.begin(),around_the_first_dig_point.end(),i)){
            continue;
        }else{
            decide_where_is_bomb.push_back(i);
        }
    }
    random_device seed_gen;
    mt19937 engine(seed_gen());
    shuffle(decide_where_is_bomb.begin(), decide_where_is_bomb.end(), engine);
    for(i = 0;i < bombnumber;++i){
        data_board_vector.at(decide_where_is_bomb.at(i)/columnnumber).at(decide_where_is_bomb.at(i)%columnnumber) = -1;
    }
    count_around_bomb(data_board_vector, linenumber, columnnumber);
    find_blank_space(display_board_vector, data_board_vector, line_dig_point, column_dig_point);
}

void count_around_bomb(vector<vector<int>> &data_board_vector, int linenumber, int columnnumber)
{
    int line,column;
    int i,j;
    for(line = 0;line < linenumber;++line){
        for(column = 0;column < columnnumber;++column){
            if(data_board_vector.at(line).at(column) == -1){
                continue;
            }else{
                for(i = -1;i <=1;++i){
                    for(j = -1;j <= 1;++j){
                        if(i == 0 && j == 0){
                            continue;
                        }else{
                            try{
                                if(data_board_vector.at(line+i).at(column+j) == -1){
                                    ++data_board_vector.at(line).at(column);
                                }
                            }catch(exception& e){
                                continue;
                            }
                        }
                    }
                }
            } 
        }
    }
}

void find_blank_space(vector<vector<char>> &display_board_vector, vector<vector<int>> &data_board_vector, int line_dig_point, int column_dig_point)
{
    queue<pair<int,int>> find_blank;
    int line,column;
    find_blank.emplace(line_dig_point,column_dig_point);
    int i,j;
    while(!find_blank.empty()){
        line = find_blank.front().first;
        column = find_blank.front().second;
        for(i = -1;i <= 1;++i){
            for(j = -1;j <= 1;++j){
                try{
                    if(display_board_vector.at(line+i).at(column+j) == 'o'){
                        if(data_board_vector.at(line+i).at(column+j) == 0){
                            display_board_vector.at(line+i).at(column+j) = ' ';
                            find_blank.emplace(line+i,column+j);
                        }else{
                            display_board_vector.at(line+i).at(column+j) = '0' + data_board_vector.at(line+i).at(column+j);
                        }
                    }else{
                        continue;
                    }
                }catch(exception &e){
                    continue;
                }
            }   
        }
        find_blank.pop();
    }
}

bool judge_game_clear(vector<vector<char>> &display_board_vector, vector<vector<int>> &data_board_vector, int linenumber, int columnnumber)
{
    int line,column;
    for(line = 0;line < linenumber;++line){
        for(column = 0;column < columnnumber;++column){
            if(data_board_vector.at(line).at(column) != -1){
                if(display_board_vector.at(line).at(column) == 'o' || display_board_vector.at(line).at(column) == 'F'){
                    return false;
                }else{
                    continue;
                }
            }else{
                continue;
            }
        }
    }
    return true;
}

int know_digit(int number)
{
    int digit = 0;
    while(number){
        number /= 10;
        ++digit;
    }
    return digit;
}