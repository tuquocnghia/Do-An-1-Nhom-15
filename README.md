# Đồ Án 1: Ma Trận và Tính Toán Khoa Học

**Môn học:** Toán Ứng Dụng và Thống Kê (MTH00051)  
**Nhóm thực hiện:** Nhóm 15  

---

## Giới thiệu đồ án
Đồ án này tập trung giải quyết các bài toán như:
* Giải hệ phương trình tuyến tính (Nghiệm duy nhất, Vô nghiệm, Vô số nghiệm).
* Tính định thức của ma trận vuông.
* Tìm ma trận nghịch đảo.
* Trích xuất hạng và cơ sở của các không gian (Cột, Dòng, Nghiệm).

*Lưu ý: Thư viện `NumPy` chỉ được sử dụng để kiểm chứng kết quả.*

---

## Thành viên nhóm và Phân công

* **21120449 - Nguyễn Văn Hậu:** Xây dựng thuật toán cốt lõi (`gaussian.py`) đưa ma trận về RREF với Partial Pivoting.
* **24120389 - Từ Quốc Nghĩa:** Xử lý biện luận nghiệm, giải hệ tam giác và xuất chuỗi công thức nghiệm tổng quát (`gaussian.py`).
* **24120263 - Đỗ Ngọc Gia Bảo:** Cài đặt các ứng dụng đại số (`determinant.py`, `inverse.py`, `rank_basis.py`).
* **24120430 - Trần Thanh Sơn:** Kiểm thử, ghi Docstring, xây dựng hệ thống Test Cases đối chiếu với NumPy, review và sửa code tạo file demo(.ipynb), đóng gói và cấu trúc thư mục đồ án.
* **23120027 - Nguyễn Hải Đăng:** Tổng hợp, viết báo cáo LaTeX.

---

## Cấu trúc thư mục


```text
Đồ án 1 (Nhóm 15)
├── README.md
├── requirements.txt
├── report/
│   └── report.pdf
└── part1/
    ├── gaussian.py         # Chứa thuật toán khử Gauss, thế ngược và biện luận nghiệm
    ├── determinant.py      # Hàm tính định thức
    ├── inverse.py          # Hàm tìm ma trận nghịch đảo
    ├── rank_basis.py       # Hàm tìm hạng và cơ sở
    └── part1_demo.ipynb    # File Jupyter Notebook để chạy Demo và Test
```