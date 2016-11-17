initial_matrices = '''
\[
A  =
\\begin{bmatrix}
%(matrice_A)s
\end{bmatrix}
\]
\\\\
\[ Initial\ set\  of\  basic\  and\  nonbasic\  indices \]
\[
\\beta= \left\{%(b_indices)s\\right\} \quad and \quad  \\nu=\left\{%(non_b_indices)s\\right\}
\]
\\\\
\[
Submatrice\ of \ A
\]

\[
B =
\\begin{bmatrix}
%(matrice_B)s
\end{bmatrix} \quad and \quad
\mathit{N} =
\\begin{bmatrix}
%(matrice_N)s
\end{bmatrix}
\]
\\\\
\[
Inital\ values\ of\ the\ basic\ variables\ are\ given\ by
\]
\[
x_B^* = b =
\\begin{bmatrix}
%(matrice_X)s
\end{bmatrix}
\]
\[
Inital\ values\ of\ the\ nonbasic\ dual variables\ are\ given\ by
\]
\[
z_\mathit{N}^*= -c_\mathit{N} =
\\begin{bmatrix}
%(matrice_Z)s
\end{bmatrix}
\]
\\\\
'''

old = '''\documentclass [12pt] {article}
\usepackage{amsmath}
\makeatletter
\\renewcommand{\@seccntformat}[1]{}
\makeatother
\usepackage{url}
\usepackage[margin=0.8in]{geometry}
\pagestyle{plain}
\\begin{document}
\section*{Initial Matrices and Vector}

\[
A  =
\\begin{bmatrix}
%(matrice_a)s
\end{bmatrix}
\]

\end{document}
'''
initial_primal_condition_true = '''\[
Since\ x_B^*\  \geq \  0, the\ initial\ solution\ is\ primal\ feasible.
\]
'''
initial_primal_condition_false = '''\[
Since\ x_B^*\  < \  0, the\ initial\ solution\ is\ primal\ feasible.
\]
'''

step1_primal_condition_true = '''\section*{Iteration No %(iteration_no)s}
\subsection{Step 1.}
\[
Since\ z_\mathit{N}^*\ has\ some\ negative\ components,\ the\ current\ solution\ is\ not\ optimal.
\]
'''

step1_primal_condition_false = '''\section*{Iteration No %(iteration_no)s}
\subsection{Step 1.}
\[
Since\ z_\mathit{N}^*\ has\ all\ nonnegative\ components,\ the\ current\ solution\ is\ optimal.
\]
'''
step2_primal = '''\subsection{Step 2.}
\[
Since\ z_\mathit{N}^*\ = \ %(negative_no)s\ and\ this\ is\ the\ most\ negative\ dual\ variables,
\]

\[
we\ see\ that\ the\ entering\ index\ is\  j\ =\ %(index)s
\]
'''

step3_primal = '''\subsection{Step 3.}
\[
\Delta x_{\mathcal B} = B^{-1} N e_j =
\\begin{bmatrix}
%(matrice_BN)s
\end{bmatrix}
\\begin{bmatrix}
%(matrice_EJ)s
\end{bmatrix}
= \\begin{bmatrix}
%(matrice_Result)s
\end{bmatrix}
\]
'''



from latex_helper import latexMatrice
from latex_helper import latexCommaSeprated
from latex_helper import latexRatios
from latex_helper import latexFraction
from latex_helper import latexVariablesMatrice

def getInitialMatrices(matrice_A, b_indices, non_b_indices, matrice_B, matrice_N, matrice_X, matrice_Z ):
    matrice_A = latexMatrice(matrice_A)
    b_indices = latexCommaSeprated(b_indices)
    non_b_indices = latexCommaSeprated(non_b_indices)
    matrice_B = latexMatrice(matrice_B)
    matrice_N = latexMatrice(matrice_N)
    matrice_X = latexMatrice(matrice_X)
    matrice_Z = latexMatrice(matrice_Z)

    text = (initial_matrices % {'matrice_A':matrice_A, 'b_indices':b_indices, 'non_b_indices':non_b_indices, 'matrice_B':matrice_B, 'matrice_N':matrice_N, 'matrice_X':matrice_X, 'matrice_Z':matrice_Z})
    return text

# return string based on whether intial condition met (i.e. true else false)
def getPrimalInitialCondition(bool_value):
    if bool_value:
        return initial_primal_condition_true
    else:
        return initial_primal_condition_false

# bool_value is false in case z_n has some negative; else true and algo stops
# iteration_no is the iteration number fo the algo
def firstStepPrimal(iteration_no,bool_value):
    text = ''
    if bool_value:
        text = (step1_primal_condition_true % {'iteration_no':iteration_no})
    else:
        text = (step1_primal_condition_false % {'iteration_no':iteration_no})
    return text


# second step of primal
# argument 1 is the most negative number
# index is the index of the most negative in the nonbasic vector
def secondStepPrimal(negative_no, index):
    text = (step2_primal % {'negative_no':negative_no,'index':index})
    return text

# returns the latex string for the third Step of primal method
def thirdStepPrimal(matrice_BN, matrice_EJ,matrice_Result):

    matrice_BN = latexMatrice(matrice_BN)
    matrice_EJ = latexMatrice(matrice_EJ)
    matrice_Result = latexMatrice(matrice_Result)

    text = (step3_primal % {'matrice_BN':matrice_BN, 'matrice_EJ':matrice_EJ,'matrice_Result':matrice_Result})
    return text

step4_primal = '''\subsection{Step 4.}
\[
t =\Bigg(
max= \left\{%(ratio)s\\right\}
\Bigg)^{-1}\ =\ %(value)s
\]
'''
def fourthStepPrimal(num_array,den_array,value):
    ratio = latexRatios(num_array,den_array)
    text = (step4_primal % {'ratio':ratio,'value':value})

    return text

step5_primal = '''\subsection{Step 5.}
\[
In\ step\ 4, \ the\ ratio\ corresponds\ to\ basic\ index\ %(index)s
\]
\[
i\ = \ %(index)s
\]
'''

def fifthStepPrimal(index):
    text = (step5_primal %{'index':index})
    return text

step6_primal='''\subsection{Step 6.}
\[
\Delta z_{\mathcal N}= -( B^{-1} N )^{T}e_i = -\
\\begin{bmatrix}
%(matrice_BN)s
\end{bmatrix}
\\begin{bmatrix}
%(matrice_EI)s
\end{bmatrix}
= \\begin{bmatrix}
%(matrice_Result)s
\end{bmatrix}
\]
'''

def sixthStepPrimal(matrice_BN,matrice_EI,matrice_Result):
    matrice_BN = latexMatrice(matrice_BN)
    matrice_EI = latexMatrice(matrice_EI)
    matrice_Result = latexMatrice(matrice_Result)

    text = (step6_primal %{'matrice_BN':matrice_BN,'matrice_EI':matrice_EI,'matrice_Result':matrice_Result})
    return text

step7_primal = '''\subsection{Step 7.}
\[
s \ =\ \\frac{z_{\mathcal N}^{*}}{ \Delta z_{\mathcal N}}\ =\ %(ratio)s\ =\ %(value)s
\]
'''
def seventhStepPrimal(value,num,den):
    text = (step7_primal %{'ratio':latexFraction(num,den),'value':value})
    return text



step8_primal = '''\subsection{Step 8.}
\[
x_{%(i_index)s}^{*}\ =\ %(i_value)s, \quad x_{\mathcal B}^{*}\ =\
\\begin{bmatrix}
%(matrice_x_old)s
\end{bmatrix}\ -%(i_value)s\
\\begin{bmatrix}
%(matrice_x_delta)s
\end{bmatrix}\ =\
\\begin{bmatrix}
%(matrice_x_new)s
\end{bmatrix}\ ,
\]
\\

\[
z_{%(j_index)s}^{*}\ =\ %(j_value)s, \quad z_{\mathcal N}^{*}\ =\
\\begin{bmatrix}
%(matrice_z_old)s
\end{bmatrix}\ -%(j_value)s\
\\begin{bmatrix}
%(matrice_z_delta)s
\end{bmatrix}\ =\
\\begin{bmatrix}
%(matrice_z_new)s
\end{bmatrix}\ ,
\]
'''

# Following parameters passed
# j,i,old_X_b,old_Z_n,delta_X_b,delta_Z_n,self.X_b,self.Z_n
def eightStepPrimal(j,i,t,s,old_X_b,old_Z_n,delta_X_b,delta_Z_n,new_X_b,new_Z_n):
    matrice_x_old = latexMatrice(old_X_b)
    matrice_z_old = latexMatrice(old_Z_n)

    matrice_x_delta = latexMatrice(delta_X_b)
    matrice_z_delta = latexMatrice(delta_Z_n)

    matrice_x_new = latexMatrice(new_X_b)
    matrice_z_new = latexMatrice(new_Z_n)

    text = (step8_primal % {'i_index':j,'j_index':i,'i_value':t,'j_value':s,'matrice_x_old':matrice_x_old,'matrice_z_old':matrice_z_old,'matrice_x_delta':matrice_x_delta,'matrice_z_delta':matrice_z_delta,'matrice_x_old':matrice_x_old,'matrice_z_old':matrice_z_old,'matrice_x_new':matrice_x_new,'matrice_z_new':matrice_z_new})
    return text

step9_primal='''\subsection{Step 9.}

\[ New\ set\  of\  basic\  and\  nonbasic\  indices \]
\[
\\beta= \left\{%(basic)s\\right\} \quad and \quad  \\nu=\left\{%(nonbasic)s\\right\}
\]

\[
Corresponding\ new\ basis\ and\ nonbasis\ submatrices\ of\ A,
\]
\[
B =
\\begin{bmatrix}
%(matrice_B)s
\end{bmatrix} \quad and \quad
\mathit{N} =
\\begin{bmatrix}
%(matrice_N)s
\end{bmatrix}
\]

\[
New\ Basic\ primal\ variables\ and\ nonbasic\ dual\ variables:
\]
\[
x_{\mathcal B}^{*}\ =\
\\begin{bmatrix}
%(matrice_var_x)s
\end{bmatrix}\ =\
\\begin{bmatrix}
%(matrice_X)s
\end{bmatrix}\
\quad
z_{\mathcal N}^{*}\ =\
\\begin{bmatrix}
%(matrice_var_z)s
\end{bmatrix}\ =\
\\begin{bmatrix}
%(matrice_Z)s
\end{bmatrix}
\]
'''

def ninthStepPrimal(basic,nonbasic,matrice_B,matrice_N,matrice_X,matrice_Z):
    matrice_var_x = latexVariablesMatrice('x',basic)
    matrice_var_z = latexVariablesMatrice('z',nonbasic)

    basic = latexCommaSeprated(basic)
    nonbasic = latexCommaSeprated(nonbasic)
    matrice_B = latexMatrice(matrice_B)
    matrice_N = latexMatrice(matrice_N)

    matrice_X = latexMatrice(matrice_X)
    matrice_Z = latexMatrice(matrice_Z)

    text = (step9_primal % {'basic':basic,'nonbasic':nonbasic,'matrice_B':matrice_B,'matrice_N':matrice_N,'matrice_X':matrice_X,'matrice_Z':matrice_Z,'matrice_var_x':matrice_var_x,'matrice_var_z':matrice_var_z})
    return text


final_doc = '''\documentclass [12pt] {article}
\usepackage{amsmath}
\makeatletter
\\renewcommand{\@seccntformat}[1]{}
\makeatother
\usepackage{url}
\usepackage[margin=0.8in]{geometry}
\pagestyle{plain}
\\begin{document}
\section*{Simplex Method Initial Matrices and Vector} %(latex_tex)s \end{document}
'''

def getWholeLatex(latex_text):
    text = (final_doc % {'latex_tex':latex_text})
    return text


objectiveFuntion ='''
\[
\zeta^{*} = %(equation)s
\]
'''
def latexObjectiveFuntion(vec_c,value):
    equation = ''
    for i in range(len(vec_c)):
        if i != len(vec_c)-1:
            equation += "%sx_{%s}^{*}\\ +" % (str(vec_c[i]),str(i+1))
        else:
            # equation += "%sx_{%s}^{*}\\ +" % (str(vec_c[i]),str(i+1))
            equation += "%sx_{%s}^{*}\\ =\\ %s" % (str(vec_c[i]),str(i+1),str(value))

    text = (objectiveFuntion %{'equation':equation})
    return text



















