from gaussian import gaussian_eliminate

def inverse(A):
    """
    [TÌM MA TRẬN NGHỊCH ĐẢO] Giải ma trận [A|I] bằng Gauss-Jordan.
    Lấy nửa bên phải của RREF làm ma trận nghịch đảo.

    Args:
        A (list of lists): Ma trận vuông cần tìm nghịch đảo.

    Returns:
        list of lists or str: Ma trận nghịch đảo (list of lists) nếu ma trận khả nghịch, 
        hoặc chuỗi lỗi ("Ma tran khong vuong" hoặc "Ma tran suy bien, khong co nghich dao") nếu không thể tìm nghịch đảo.
    """
    m = len(A)
    n = len(A[0])
    if m != n:
        return "Ma tran khong vuong"

    I = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    RREF_matrix, _, _, det_multiplier = gaussian_eliminate(A, I)

    if abs(det_multiplier) < 1e-9:
        return "Ma tran suy bien, khong co nghich dao"

    inv_matrix = []
    for r in range(m):
        inv_matrix.append(RREF_matrix[r][n:])

    return inv_matrix