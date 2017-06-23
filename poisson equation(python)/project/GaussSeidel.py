def Gauss(matrix,value):
    """一次高斯迭代更新格点上的值"""
    value = value.copy()
    for i in range(len(matrix)):
        result = 0
        for j in matrix[i]:
            if i != j:
                result += matrix[i][j] * value[j]
        value[i] = (self.get_func_value(i) - result) / self.matrix[i][i]
    return value
