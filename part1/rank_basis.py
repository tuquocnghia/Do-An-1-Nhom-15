from gaussian import gaussian_eliminate

def rank_and_basis(A):
    """
    [HẠNG VÀ CƠ SỞ] Tính hạng và trích xuất cơ sở Không gian Cột, Dòng, Nghiệm.

    Args:
        A (list of lists): Ma trận cần phân tích.

    Returns:
        dict: Từ điển chứa các khóa sau:
            - 'rank' (int): Hạng của ma trận.
            - 'column_space_basis' (list of lists): Cơ sở của không gian cột (danh sách các vector cột pivot).
            - 'row_space_basis' (list of lists): Cơ sở của không gian dòng (danh sách các hàng không zero trong RREF).
            - 'null_space_basis' (list of lists): Cơ sở của không gian nghiệm (danh sách các vector nghiệm cơ bản).
    """
    m = len(A)
    n = len(A[0])
    epsilon = 1e-9

    b_dummy = [0.0] * m
    RREF_matrix, _, _, _ = gaussian_eliminate(A, b_dummy)

    pivot_cols = []
    r_idx = 0
    for c_idx in range(n):
        if r_idx >= m:
            break
        if abs(RREF_matrix[r_idx][c_idx]) > epsilon:
            pivot_cols.append(c_idx)
            r_idx += 1

    rank = len(pivot_cols)
    column_space = [[A[r][c] for r in range(m)] for c in pivot_cols]

    row_space = []
    for r in range(m):
        if any(abs(val) > epsilon for val in RREF_matrix[r][:n]):
            row_space.append(RREF_matrix[r][:n])

    free_cols = [c for c in range(n) if c not in pivot_cols]
    null_space = []
    for fc in free_cols:
        vec = [0.0] * n
        vec[fc] = 1.0
        for idx, pc in enumerate(pivot_cols):
            vec[pc] = -RREF_matrix[idx][fc]
        null_space.append(vec)

    return {
        'rank': rank,
        'column_space_basis': column_space,
        'row_space_basis': row_space,
        'null_space_basis': null_space
    }