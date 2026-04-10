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
    m = len(c)
    solution = [0.0] * m
    for i in range(m - 1, -1, -1): # Duyệt từ hàng dưới cùng lên
        if abs(U[i][i]) < EPSILON: # Ứng với mỗi hàng, kiểm tra phần tử thuộc đường chéo chính có = 0 hay không. Nếu = 0 thì kết luận phương trình không có nghiệm duy nhất
            raise ValueError("He phuong trinh khong co nghiem duy nhat!")
        total = 0 # Ngược lai nếu khác 0, ta cần tính tổng các hạng tử nằm bên phải vị trí đang xét
        for j in range(i + 1, m):
            total += U[i][j] * solution[j]
        solution[i] = (c[i] - total) / U[i][i] # chuyển vế để suy ra nghiệm ở cột hiện tại
    return solution 


def solve_system(RREF_matrix):
    if not RREF_matrix or not RREF_matrix[0]: # Nếu ma trận không tồn tại thì trả về rỗng
        return "Ma tran rong"
    res_str = "" # Khởi tạo chuỗi kết quả 
    rows = len(RREF_matrix) # Số dòng ma trận
    cols = len(RREF_matrix[0]) # Số cột ma trận
    n = cols - 1

    for i in range(rows): # Duyệt từng hàng
        isAllZero = all(abs(RREF_matrix[i][j]) < EPSILON for j in range(n)) # Kiểm tra tất cả hệ số của ẩn ở hàng đó có = 0 hết hay không
        if (isAllZero and abs(RREF_matrix[i][n]) > EPSILON): # Nếu tất cả hệ số của ẩn ở hàng đó = 0 và hệ số tự do khác 0 thì kết luận vô nghiệm
            return "Phuong trinh vo nghiem"
    

    pivots = {} # Khởi tạo danh sách pivot (phần tử cơ sở)
    for i in range(rows): # Duyệt hàng
        for j in range(n): # Duyệt cột
            if abs(RREF_matrix[i][j]) > EPSILON: # Vị trí đầu tiên ở mỗi hàng khác 0 sẽ là pivot
                pivots[i] = j
                break # Tìm sang hàng tiếp theo khi đã tim được pivot cho hàng trước đó
    
    if len(pivots) == n: # Tổng số phần tử pivot = tổng số nghiệm thì hệ có nghiệm duy nhất
        solution = [0.0] * n # Tạo danh sách lưu các ẩn
        for r, c in pivots.items(): # Duyệt từng pivot, ứng với mỗi pivot thì ẩn của nó = hệ số tự do của chính hàng đó
            solution[c] = RREF_matrix[r][n] # lưu ẩn

        # Nối chuỗi định dạng (x1, x2, ..., xn) = (...)
        var_names = ", ".join([f"x_{i+1}" for i in range(n)])
        sol_values = ", ".join([f"{val:.4g}" for val in solution])
        
        res_str = f"({var_names}) = ({sol_values})"
        return f"He co nghiem duy nhat: {res_str}"
    


    free_vars = [j for j in range(n) if not j in pivots.values()] # Tìm các cột không chứa pivot, đây chính là các biến tự do
    free_var_symbols = {value: f"t_{idx + 1}" for idx, value in enumerate(free_vars)} # Đặt tên cho các biến tự do thành các tham số t_1, t_2,... để dễ biểu diễn nghiệm tổng quát

    for j in range(n): # Duyệt qua từng biến x_1, x_2, ... x_n để xuất công thức nghiệm
        if not j in pivots.values(): # Nếu cột j là biến tự do (không có pivot)
            res_str += f"x_{j+1} = {free_var_symbols[j]}\n" # Thì gán trực tiếp biến đó bằng tham số tương ứng
        else:
            # Ngược lại, nếu cột j là biến cơ sở (có chứa pivot)
            # Tìm dòng đang chứa phần tử pivot của biến x_j hiện tại
            currentRow = [r for r, c in pivots.items() if c == j][0]
            constant = RREF_matrix[currentRow][n] # Lấy hệ số tự do (nằm ở cột cuối cùng bên phải) của dòng đó
            eq_parts = [] # Khởi tạo mảng lưu các thành phần của vế phải
            if abs(constant) > EPSILON: # Nếu hệ số tự do khác 0
                eq_parts.append(f"{constant:.4g}") # Đưa hệ số tự do vào vế phải
            for free_var in free_vars: # Duyệt qua các biến tự do nằm trên cùng dòng này để chuyển vế
                coef = RREF_matrix[currentRow][free_var] # Lấy hệ số của biến tự do
                coef = -coef # Đổi dấu hệ số vì ta đang chuyển vế nó từ trái (chứa ẩn) sang phải (chứa kết quả)
                if abs(coef) > EPSILON: # Nếu hệ số sau khi chuyển vế khác 0 thì mới đưa vào công thức
                    sign = "+ " if coef >= EPSILON else "- " # Xác định dấu của hạng tử
                    coef = abs(coef) # Lấy trị tuyệt đối để dễ in format ghép với dấu ở trên
                    if abs(coef - 1) < EPSILON: # Nếu hệ số là 1
                        part = f"{sign} {free_var_symbols[free_var]}" # Tránh in số 1 cho đẹp
                    else:
                        part = f"{sign}{coef:.4g}{free_var_symbols[free_var]}" # Còn khác 1 thì in bình thường: dấu + hệ số + tên tham số
                    if not eq_parts and sign == "+ ": # Nếu đây là phần tử đứng đầu tiên của vế phải và mang dấu dương
                        part = part.strip("+ ") # Thì xóa dấu cộng đi
                    eq_parts.append(part.strip()) # Thêm hạng tử đã format gọn gàng vào mảng
            if not eq_parts:
                eq_parts.append("0") # Thì vế phải mặc định là 0
            res_str += f"x_{j+1} = {' '.join(eq_parts)}\n" # Nối các hạng tử trong mảng bằng dấu cách để tạo thành phương trình hoàn chỉnh cho x_j

    return res_str
