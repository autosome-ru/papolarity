m = [6,7,5,8,6,3,0,0,0,0,0,0,0,0]
#создали массив и присвоили ему имя m
total_m = 0
#переменная, которая равна 0. Потом будет переменной, которая в знаменатель
numerator = 0
#переменная для числителя в формуле
for i in range(len(m)):
#for - цикл; range генерирует значения от 0 (например) заданное количество значений (10, в данном случае); len(m) - считает длину массива
    total_m = total_m + m[i]
    #новое значение total_m равно предудыщее значение + i-тое значение массива m
    numerator += m[i] * i
    #
polarity_score = numerator / total_m
#вычислили "центр тяжести"
print(polarity_score)
normalized_score = (polarity_score * 2 / (len(m) - 1)) - 1
print(normalized_score)


def calculate_polarity_score(m):
#define; m - временное рандомное имя
    total_m = 0
    #переменная, которая равна 0. Потом будет переменной, которая в знаменатель
    numerator = 0
    #переменная для числителя в формуле
    for i in range(len(m)):
    #for - цикл; range генерирует значения от 0 (например) заданное количество значений (10, в данном случае); len(m) - считает длину массива
        total_m = total_m + m[i]
        #новое значение total_m равно предудыщее значение + i-тое значение массива m
        numerator += m[i] * i
        #
    polarity_score = numerator / total_m
    #вычислили "центр тяжести"
    normalized_score = (polarity_score * 2 / (len(m) - 1)) - 1
    return normalized_score
    #возвращает в программу полученное число, которым теперь она сможет пользоваться

pol = calculate_polarity_score([1,10,12,8,3,2,2,0,0])
print(pol)
#calculate_polarity_score(m)