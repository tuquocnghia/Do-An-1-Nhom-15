EPSILON = 1e-9

def gaussian_eliminate(A, b):
    """
    Đưa ma trận [A|b] về RREF bằng Gauss-Jordan (có Partial Pivoting).

    Args:
        A (list of lists): Ma trận hệ số.
        b (list/list of lists): Vector hệ số tự do hoặc ma trận ghép (VD: ma trận I).

    Returns:
        Trả về tuple (RREF_matrix, x, swap_count, det_multiplier):
        - RREF_matrix: Ma trận RREF.
        - x: Nghiệm hệ phương trình (list) hoặc chuỗi báo lỗi ("Vo nghiem"...).
        - swap_count: Số lần hoán đổi dòng.
        - det_multiplier: Tích các pivot (dùng để tính định thức).
    """
    m = len(A)
    n = len(A[0])
    epsilon = 1e-9

    M = []
    if isinstance(b[0], list):
        b_cols = len(b[0])
    else:
        b_cols = 1
    is_b_matrix = (b_cols > 1)

    for r in range(m):
        row = [float(val) for val in A[r]]
        if isinstance(b[r], list):
            row.extend([float(val) for val in b[r]])
        else:
            row.append(float(b[r]))
        M.append(row)

    swap_count = 0
    det_multiplier = 1.0
    i = 0

    for k in range(n):
        if i >= m:
            break

        max_row = i
        max_val = abs(M[i][k])
        for r in range(i + 1, m):
            if abs(M[r][k]) > max_val:
                max_val = abs(M[r][k])
                max_row = r

        if max_val < epsilon:
            print(f"Khong co pivot tai cot {k}")
            det_multiplier = 0.0
            continue

        if max_row != i:
            M[i], M[max_row] = M[max_row], M[i]
            swap_count += 1

        pivot = M[i][k]
        det_multiplier *= pivot

        for c in range(k, len(M[0])):
            M[i][c] /= pivot

        for r in range(m):
            if r != i:
                factor = M[r][k]
                for c in range(k, len(M[0])):
                    M[r][c] -= factor * M[i][c]
        i += 1

    RREF_matrix = M
    x = "Khong xac dinh"
    if not is_b_matrix:
        is_inconsistent = False
        for r in range(m):
            all_zero_left = all(abs(M[r][c]) < epsilon for c in range(n))
            if all_zero_left and abs(M[r][n]) >= epsilon:
                is_inconsistent = True
                break

        if is_inconsistent:
            x = "Vo nghiem"
        elif i < n: 
            x = "Vo so nghiem"
        else:
            x = [M[r][n] for r in range(n)]

    return (RREF_matrix, x, swap_count, det_multiplier)

def back_substitution(U, c):
    """
    Thực hiện phép thế ngược trên ma trận tam giác trên U để giải hệ Ux = c.

    Args:
        U (list of lists): Ma trận tam giác trên.
        c (list): Vector hệ số tự do.

    Returns:
        list: Nghiệm của hệ phương trình.

    Raises:
        ValueError: Nếu ma trận U không khả nghịch (pivot bằng 0).
    """
    m = len(c)
    solution = [0.0] * m
    for i in range(m - 1, -1, -1):
        if abs(U[i][i]) < EPSILON:
            raise ValueError("He phuong trinh khong co nghiem duy nhat!")
        total = 0
        for j in range(i + 1, m):
            total += U[i][j] * solution[j]
        solution[i] = (c[i] - total) / U[i][i]
    return solution

def solve_system(RREF_matrix):
    """
    Phân tích ma trận RREF để xác định nghiệm của hệ phương trình tuyến tính.

    Args:
        RREF_matrix (list of lists): Ma trận đã được đưa về dạng RREF (Reduced Row Echelon Form).

    Returns:
        str: Chuỗi mô tả nghiệm hệ phương trình. Có thể là:
             - "Ma tran rong" nếu ma trận rỗng.
             - "Phuong trinh vo nghiem" nếu hệ vô nghiệm.
             - "He co nghiem duy nhat: (x_1, x_2, ...) = (val1, val2, ...)" nếu có nghiệm duy nhất.
             - Chuỗi biểu diễn nghiệm tổng quát với biến tự do nếu có vô số nghiệm.
    """
    
    if not RREF_matrix or not RREF_matrix[0]:
        return "Ma tran rong"

    res_str = ""
    rows = len(RREF_matrix)
    cols = len(RREF_matrix[0])
    n = cols - 1

    for i in range(rows):
        isAllZero = all(abs(RREF_matrix[i][j]) < EPSILON for j in range(n))
        if (isAllZero and abs(RREF_matrix[i][n]) > EPSILON):
            return "Phuong trinh vo nghiem"

    pivots = {}
    for i in range(rows):
        for j in range(n):
            if abs(RREF_matrix[i][j]) > EPSILON:
                pivots[i] = j
                break

    if len(pivots) == n:
        solution = [0.0] * n
        for r, c in pivots.items():
            solution[c] = RREF_matrix[r][n]
        var_names = ", ".join([f"x_{i + 1}" for i in range(n)])
        sol_values = ", ".join([f"{val:.4g}" for val in solution])
        return f"He co nghiem duy nhat: ({var_names}) = ({sol_values})"

    free_vars = [j for j in range(n) if j not in pivots.values()]
    free_var_symbols = {val: f"t_{idx + 1}" for idx, val in enumerate(free_vars)}

    for j in range(n):
        if j not in pivots.values():
            res_str += f"x_{j + 1} = {free_var_symbols[j]}\n"
        else:
            currentRow = [r for r, c in pivots.items() if c == j][0]
            constant = RREF_matrix[currentRow][n]
            eq_parts = []
            if abs(constant) > EPSILON:
                eq_parts.append(f"{constant:.4g}")
            for free_var in free_vars:
                coef = -RREF_matrix[currentRow][free_var]
                if abs(coef) > EPSILON:
                    sign = "+ " if coef > 0 else "- "
                    coef = abs(coef)
                    part = f"{sign}{free_var_symbols[free_var]}" if abs(coef - 1) < EPSILON else f"{sign}{coef:.4g}{free_var_symbols[free_var]}"
                    if not eq_parts and sign == "+ ":
                        part = part.strip("+ ")
                    eq_parts.append(part.strip())
            if not eq_parts:
                eq_parts.append("0")
            res_str += f"x_{j + 1} = {' '.join(eq_parts)}\n"

    return res_str