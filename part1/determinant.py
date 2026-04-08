from gaussian import gaussian_eliminate

def determinant(A):
    """
    [TÍNH ĐỊNH THỨC] Gọi hàm Gauss để lấy tích các pivot và số lần hoán đổi.
    Công thức: det(A) = (-1)^swap_count * det_multiplier

    Args:
        A (list of lists): Ma trận vuông cần tính định thức.

    Returns:
        float or str: Giá trị định thức (float) nếu ma trận vuông và khả nghịch, hoặc chuỗi "Ma tran khong vuong, khong tinh duoc dinh thuc" nếu ma trận không vuông.
    """
    m = len(A)
    n = len(A[0])
    if m != n:
        return "Ma tran khong vuong, khong tinh duoc dinh thuc"

    b_dummy = [0.0] * m
    _, _, swap_count, det_multiplier = gaussian_eliminate(A, b_dummy)

    det = ((-1) ** swap_count) * det_multiplier
    return 0.0 if abs(det) < 1e-9 else det