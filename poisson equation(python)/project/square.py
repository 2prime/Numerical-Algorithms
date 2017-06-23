'''
方阵模块
'''
__version__ = 1.3


def listcopy(lst):
    result = []
    for i in lst:
        result.append(i.copy())
    return result


def Gauss_equation(Gauss_matrix, args):
    answer = []
    for RHS in args:
        RHS = RHS.copy()
        assert len(RHS) == Gauss_matrix.order, "matrix的阶数应该等于RHS的长度"
        for j in range(Gauss_matrix.order -1, 0, -1):
            RHS[j] = RHS[j] / Gauss_matrix.content[j][j]
            for i in range(j):
                RHS[i] -= RHS[j] * Gauss_matrix.content[i][j]
        RHS[0] /= Gauss_matrix.content[0][0]
        answer.append(RHS)
    if len(answer) == 1:
        return answer[0]
    return answer



class square_matrix():
    def __init__(self,content):
        self.content = listcopy(content)
        self.order = len(content[0])
        self.col_num = self.order
        self.row_num = self.order

    def ele(self,x,y):
        assert 1<=x<=self.order and 1<=y<=self.order
        return self.content[x-1][y-1]

    def __add__(self, other):
        assert self.order == other.order,"两个不相同阶的矩阵不能相加"
        result = [[0 for i in range(self.order)] for i in range(self.order)]
        for i in range(self.order):
            for j in range(self.order):
                result[i][j] = self.content[i][j] + other.content[i][j]
        return square_matrix(result)

    def conj(self):
        return square_matrix(list(map(list, list(zip(*self.content)))))

        result = [[0 for i in range(self.order)] for i in range(self.order)]
        for i in range(self.order):
            for j in range(self.order):
                result[i][j] = self.content[j][i]
        return square_matrix(result)


    def change(self,m,n):
        for i in range(self.order):
            self.content[m-1][i],self.content[n-1][i] = self.content[n-1][i],self.content[m-1][i]
        # for i in range(self.order):
            # self.content[i][m - 1], self.content[i][n - 1] = self.content[i][n - 1], self.content[i][m - 1]

    def Gauss(self,args=[]):
        flag = not args == []
        Gauss_result = listcopy(self.content)

        change_times = 0

        def lst_change(lst,m, n):
            for i in range(self.order):
                lst[m - 1][i], lst[n - 1][i] = lst[n - 1][i], lst[m - 1][i]

        for times in range(1,self.order):
            A = square_matrix(Gauss_result)
            # print(A)
            correct_cof = {}
            max_position = max([i for i in range(times,self.order+1)],key=lambda x:abs(A.ele(x,times)))
            # print(max_position)
            A.change(max_position,times)
            lst_change(Gauss_result,max_position,times)
            if max_position != times:
                change_times += 1
            # print(A)
            # print(result)
            if flag:
                for result in args:
                    result[max_position-1],result[times-1] = result[times-1],result[max_position-1]
            # print(result)
            for i in range(times+1,self.order+1):
                correct_cof[i] = A.ele(i,times)/A.ele(times,times)
            for i in range(times+1,self.order+1):
                for j in range(times,self.order+1):
                    Gauss_result[i - 1][j - 1] = A.ele(i, j) - correct_cof[i] * A.ele(times, j)
                if flag:
                    for result in args:
                        result[i-1] -= correct_cof[i] * result[times-1]
            # print(result)
            # print(A)
            # print("======")
        if flag:
            # print("end")
            return [square_matrix(Gauss_result)]+list(args)
        return square_matrix(Gauss_result)

    def det(self):
        Gauss_result = listcopy(self.content)
        change_times = 0

        def lst_change(lst, m, n):
            for i in range(self.order):
                lst[m - 1][i], lst[n - 1][i] = lst[n - 1][i], lst[m - 1][i]

        for times in range(1, self.order):
            A = square_matrix(Gauss_result)
            # print(A)
            correct_cof = {}
            max_position = max([i for i in range(times, self.order + 1)], key=lambda x: abs(A.ele(x, times)))
            # print(max_position)
            if not A.ele(max_position,times):
                return 0
            A.change(max_position, times)
            lst_change(Gauss_result, max_position, times)
            if max_position != times:
                change_times += 1
            # print(A)
            # print(result)
            # print(result)
            for i in range(times + 1, self.order + 1):
                correct_cof[i] = A.ele(i, times) / A.ele(times, times)
            for i in range(times + 1, self.order + 1):
                for j in range(times, self.order + 1):
                    Gauss_result[i - 1][j - 1] = A.ele(i, j) - correct_cof[i] * A.ele(times, j)
        det = 1
        for i in range(self.order):
            det *= Gauss_result[i][i]
        return det*((-1)**change_times)


    # def conj(self):
    #     return square_matrix(list(map(list, list(zip(*self.content)))))

    #    result = [[0 for i in range(self.order)] for i in range(self.order)]
    #    for i in range(self.order):
    #        for j in range(self.order):
    #            result[i][j] = self.content[j][i]
    #    return square_matrix(result)

    def inv(self):
        answer = []
        for i in range(self.order):
            need_to_cal = [0 for i in range(self.order)]
            need_to_cal[i] = 1
            answer.append(need_to_cal.copy())

        return square_matrix(self.result(answer)).conj()


    def result(self,args):
        Gauss = self.Gauss(args)
        return Gauss_equation(Gauss[0],Gauss[1:])

    def __str__(self):
        string = ""
        for i in range(self.order):
            for j in range(self.order):
                string += ("%5.10f")%(self.content[i][j])+"  "
            string += "\n"
        return string



if __name__ == "__main__":
    import random,time
    A = [[6 if i==j else 8 if i-j ==1 else 1 if i-j==-1 else 0 for j in range(84)]for i in range(84)]
    b = [7]+[15]*82+[14]
    num = 1
    #t = time.time()
    #print(square_matrix(A).det())
    #print(time.time() - t)
    #print(time.time()-t)
    #print(square_matrix(A).inv())
    if True:
        t = time.time()
        for i in square_matrix(A).result([b]):
            print("%.5f"%(i),end="  ")
            if num%10 ==0:
                print("\n")
            num += 1
        #print(time.time()-t)
    if False:
        print("\n=======")
        print(square_matrix([[1,4,7],[2,5,8],[3,6,10]]).inv())
        print("\n=======")
    if False:
        A = [[10 if i==j else 1 if i-j ==1 else 1 if i-j==-1 else 0 for j in range(100)]for i in range(100)]
        b = [[random.random()for i in range(100)]]
        num =1
        t  = time.time()
        for i in square_matrix(A).result(b):
            print("%.5f"%(i), end="  ")
            if num % 10 == 0:
                print("\n")
            num += 1
        print("\n=======")
        print(time.time()-t)
    A=[[1,2,3],[4,5,6],[7,8,9]]
    #print(square_matrix(A).det())