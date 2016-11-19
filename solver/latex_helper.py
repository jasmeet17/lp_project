import numpy as np
import json

# returns the matrice suitalble for Latex format
def latexMatrice(numpy_array):
    text = ''
    if len(numpy_array.shape)==1:
        for i in range(numpy_array.shape[0]):
            text += str(numpy_array[i]) + " \\\\ "
    else:
        for i in range(numpy_array.shape[0]):
            for j in range(numpy_array.shape[1]):
                if j<numpy_array.shape[1]-1:
                    text += str(numpy_array[i][j]) + " & "
                else:
                    text += str(numpy_array[i][j]) + " \\\\ "

    return text

# return comma seprated array example [3,4,5] -> 3,4,5
def latexCommaSeprated(numpy_array):
    text = json.dumps(numpy_array,',')
    return text[1:-1]

# return output of the form # \frac{1}{1},\frac{2}{0},\frac{0}{5}
def latexRatios(num_array,den_array):
    text = ''
    for i in range(len(num_array)):
        if i < len(num_array) - 1:
            text += '\\frac{%s}{%s},' % (num_array[i],den_array[i])
        else:
            text += '\\frac{%s}{%s}' % (num_array[i],den_array[i])

    return text

# return output of the form \frac{-4}{-1}
def latexFraction(num,den):
    text = '\\frac{%s}{%s}' % (num,den)
    return text


# c contains teh letter of variable
# arr contains teh subscipts
def latexVariablesMatrice(c , arr):
    text = ''
    for i in range(len(arr)):
        text += '%s_{%s}^{*} \\\\' %(c,arr[i])
    return text
