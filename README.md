# Compression - Lượng tử hóa đều và nén giãn (μ-law, A-law) cho tín hiệu tiếng nói

Mô phỏng trên MATLAB so sánh **lượng tử hóa đều** với **lượng tử hóa không đều thông qua nén giãn (companding)** theo hai chuẩn viễn thông: **luật μ** (μ = 255, Bắc Mỹ và Nhật Bản) và **luật A** (A = 87.6, Châu Âu). Tín hiệu thử nghiệm là một đoạn tiếng nói nam; hiệu năng được đánh giá bằng **phương sai sai số lượng tử** và **tỉ số tín hiệu trên nhiễu (SNR)** tính bằng phương pháp số.

> Đồ án môn **Truyền thông số** - HK2, năm học 2025-2026
> Khoa Điện tử - Viễn thông, Trường Đại học Khoa học Tự nhiên, ĐHQG TP.HCM
> Lớp **23DTV_CLC1** - Nhóm 8 - Giảng viên hướng dẫn: TS Đặng Lê Khoa - TP.HCM, 02/2026

## Thành viên

| STT | Họ và tên | Ghi chú |
|---|---|---|
| 1 | Lý Anh Thư | |
| 2 | Hà Vĩ Khang | |
| 3 | Đỗ Trường Tín | |
| 4 | Mai Lý Khải Triều | Nhóm trưởng |

## Mục lục

- [1. Bài toán](#1-bài-toán)
- [2. Cơ sở lý thuyết](#2-cơ-sở-lý-thuyết)
- [3. Thông số mô phỏng](#3-thông-số-mô-phỏng)
- [4. Quy trình mô phỏng](#4-quy-trình-mô-phỏng)
- [5. Cách chạy](#5-cách-chạy)
- [6. Kết quả](#6-kết-quả)
- [7. Nhận xét và lưu ý về kết quả](#7-nhận-xét-và-lưu-ý-về-kết-quả)
- [8. Hướng phát triển](#8-hướng-phát-triển)
- [9. Cấu trúc thư mục đề xuất](#9-cấu-trúc-thư-mục-đề-xuất)
- [10. Tài liệu tham khảo](#10-tài-liệu-tham-khảo)

## 1. Bài toán

Trong các hệ thống truyền thông số, lượng tử hóa đều tín hiệu giọng nói thường cho SNR thấp. Nguyên nhân là đặc tính thống kê của tiếng nói: các mẫu có biên độ nhỏ xuất hiện với tần suất cao hơn nhiều so với các mẫu biên độ lớn, trong khi lượng tử hóa đều dùng cùng một bước nhảy `q` cho toàn dải nên sai số tương đối với tín hiệu nhỏ rất lớn.

Kỹ thuật **nén giãn** nén tín hiệu theo đặc tuyến logarit trước khi lượng tử hóa rồi giãn ngược lại ở phía thu, nhằm phân bổ nhiều mức lượng tử hơn cho vùng biên độ nhỏ. Đồ án lập trình và so sánh ba hệ thống với cùng số mức lượng tử L = 16:

1. Lượng tử hóa đều.
2. Nén giãn theo luật μ.
3. Nén giãn theo luật A.

## 2. Cơ sở lý thuyết

### 2.1. Định dạng tín hiệu âm thanh (PCM)

Quá trình PCM gồm ba bước: **lấy mẫu** (tần số lấy mẫu `fs >= 2fm` theo định lý Nyquist), **lượng tử hóa** (làm tròn mẫu về mức lượng tử gần nhất) và **mã hóa nhị phân**.

| Đại lượng | Công thức |
|---|---|
| Số bit mã hóa | `l = log2(L)` |
| Tốc độ bit PCM | `Rb = fs * n` (n là số bit lượng tử) |
| Tần số chồng phổ | `f_alias = abs(f - k*fs)` |

Với `fs = 4000 Hz` và `L = 16` (4 bit mỗi mẫu), tốc độ bit là 16 kbit/s.

### 2.2. Lượng tử hóa đều

Các mức lượng tử cách đều nhau một bước nhảy `q`. Sai số lượng tử có hai loại:

- **Sai số tuyến tính:** tín hiệu nằm trong dải động của bộ lượng tử, sai số giới hạn trong `±q/2`.
- **Sai số bão hòa:** tín hiệu vượt dải động, bị cắt ngưỡng, sai số thường lớn hơn nhiều.

Phương sai nhiễu lượng tử gồm hai thành phần: `σq² = σ_lin² + σ_sat²`.

### 2.3. Lượng tử hóa không đều và nén giãn

**Luật μ (μ = 255):**

```
y = y_max * ln(1 + μ*|x|/x_max) / ln(1 + μ) * sgn(x)
```

**Luật A (A = 87.6):**

```
0 < |x|/x_max <= 1/A :  y = y_max * (A*|x|/x_max) / (1 + ln A) * sgn(x)
1/A < |x|/x_max < 1  :  y = y_max * (1 + ln(A*|x|/x_max)) / (1 + ln A) * sgn(x)
```

Phía thu dùng hàm ngược (giãn) của từng luật để khôi phục tín hiệu.

### 2.4. Chỉ tiêu đánh giá

- **Phương sai sai số lượng tử** `σq²`: trung bình bình phương của sai số giữa tín hiệu gốc và tín hiệu khôi phục, càng nhỏ càng tốt.
- **SNR** = `P_signal / P_noise`, với nhiễu chính là nhiễu lượng tử, càng lớn càng tốt.

## 3. Thông số mô phỏng

| Tham số | Giá trị |
|---|---|
| Tín hiệu vào | `MaleSpeech-16-4-mono-20secs.wav` (tiếng nói nam, mono, 20 giây) |
| Tần số lấy mẫu `Fs` | 4000 Hz |
| Đoạn xử lý | 1.5 giây đầu (`t = 0:1/Fs:1.5`, 6001 mẫu) |
| Số mức lượng tử `L` | 16 |
| Điện áp đỉnh `V_p` | 0.5625 V |
| Bước nhảy `q` | `V_p/(L-1) = 0.0375` |
| Chỉ số mức | từ -8 đến 7 (`min_index = -L/2`, `max_index = L/2 - 1`) |
| `x_max`, `y_max` | `V_p` |
| μ, A | 255, 87.6 |
| Công suất tín hiệu `pow_sig` | 0.0106 |

## 4. Quy trình mô phỏng

| Bước | Nội dung |
|---|---|
| 1 | Nạp file âm thanh, lấy 1.5 giây đầu làm tín hiệu `mSpeech` |
| 2 | Lượng tử hóa đều `mSpeech` với L = 16, thu được `s_q2` |
| 3 | Vẽ tín hiệu gốc và tín hiệu lượng tử hóa đều |
| 4 | Tính phương sai sai số `var_sq2` và SNR `snr_sq2` của lượng tử hóa đều |
| 5 | Nén `mSpeech` bằng luật μ và luật A, thu được `s_c5_mu`, `s_c5_A` |
| 6 | Lượng tử hóa tín hiệu đã nén với cùng số mức và bước nhảy như bước 2 (`s_q6_mu`, `s_q6_A`) |
| 7 | Giãn tín hiệu đã lượng tử hóa để khôi phục (`s_e7_mu`, `s_e7_A`) |
| 8 | Vẽ các tín hiệu nén, lượng tử hóa sau nén và giãn khôi phục cùng đồ thị với tín hiệu gốc (phóng to đoạn 0.5275 s đến 0.59 s) |
| 9 | Tính phương sai và SNR của hai hệ thống nén giãn, đối chiếu với lượng tử hóa đều |

## 5. Cách chạy

**Yêu cầu:** MATLAB (chỉ dùng các hàm cơ bản: `audioread`, `round`, `log`, `exp`, `sign`, `plot`).

1. Đặt file `MaleSpeech-16-4-mono-20secs.wav` cùng thư mục với script.
2. Chép mã MATLAB trong mục *3.2. Chương trình Matlab* của báo cáo (dòng 001-101) vào file `companding_simulation.m`.
3. Chạy script. Chương trình vẽ đồ thị so sánh và in kết quả ra Command Window:

```
KET QUA
1. Luong tu hoa deu (Uniform) - B4: SNR = 30.2781
2. Companding (u-Law) - B9: SNR = 1.5499
3. Companding (A-Law) - B9: SNR = 1.6207
```

Các biến chính trong Workspace: `var_sq2`, `var_mu`, `var_A` (phương sai) và `snr_sq2`, `snr_mu`, `snr_A` (SNR, dạng tỉ số tuyến tính).

## 6. Kết quả

| Hệ thống | Phương sai sai số | SNR (tỉ số) | SNR (dB) |
|---|---|---|---|
| Lượng tử hóa đều | 3.4921e-04 | 30.2781 | 14.81 |
| Nén giãn luật μ (μ = 255) | 0.0068 | 1.5499 | 1.90 |
| Nén giãn luật A (A = 87.6) | 0.0065 | 1.6207 | 2.10 |

Nếu so sánh hai luật nén giãn: luật A cho phương sai nhỏ hơn và SNR cao hơn luật μ một chút. Báo cáo giải thích rằng luật A có đoạn tuyến tính ở vùng biên độ nhỏ (`0 < |x|/x_max <= 1/A`), nên phép giãn không dùng hàm mũ ở vùng quanh gốc tọa độ, là vùng chứa nhiều mẫu tiếng nói nhất, do đó kìm được sự bùng nổ của sai số lượng tử.

## 7. Nhận xét và lưu ý về kết quả

**Lưu ý quan trọng:** trong lần chạy này, SNR của hai hệ thống nén giãn (khoảng 1.5-1.6) **thấp hơn nhiều** so với lượng tử hóa đều (30.28), nên các con số chưa cho thấy companding cải thiện SNR trung bình so với lượng tử hóa đều.

Nguyên nhân nhiều khả năng nằm ở dải động của bộ lượng tử trong mã mô phỏng:

- Với `q = V_p/(L-1) = 0.0375` và chỉ số mức từ -8 đến 7, bộ lượng tử chỉ biểu diễn được các giá trị từ `-8q = -0.3` đến `7q = 0.2625`, tức là chưa bằng nửa dải `±V_p = ±0.5625`.
- Tín hiệu sau khi nén có biên độ trải rộng đến `±V_p`, nên phần lớn mẫu bị **cắt ngưỡng (bão hòa)** khi lượng tử hóa. Các điểm lượng tử hóa sau nén trên đồ thị của báo cáo dừng ở khoảng +0.26 và -0.30, khớp với hai giá trị này.
- Sai số bão hòa lớn hơn nhiều sai số tuyến tính, làm phương sai của hai hệ thống nén giãn tăng mạnh (0.0065-0.0068, so với 3.49e-04 của lượng tử hóa đều).

Vì vậy kết quả hiện tại phản ánh ảnh hưởng của bão hòa nhiều hơn là hiệu quả thực của luật μ và luật A. Để so sánh công bằng, cần chọn `q` sao cho L mức phủ trọn dải `±V_p` (ví dụ `q = 2*V_p/L`) và dùng cùng dải động cho cả ba hệ thống, rồi chạy lại và cập nhật bảng kết quả ở mục 6.

## 8. Hướng phát triển

- Chỉnh dải động bộ lượng tử như mục 7 rồi đối chiếu lại ba hệ thống.
- Báo cáo SNR theo đơn vị dB và vẽ **SNR theo biên độ tín hiệu vào**, để thấy rõ ưu điểm của companding với tín hiệu nhỏ.
- Khảo sát ảnh hưởng của số mức L (số bit 4 đến 8) và của tham số μ, A.
- Thử nhiều file tiếng nói khác (nam, nữ, nhiều tốc độ nói).
- Đóng gói chương trình thành các hàm (`compress`, `expand`, `quantize`) để dễ tái sử dụng và kiểm thử.

