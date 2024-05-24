function [G,x,y,z,alpha] = generate_G_near_field_channel(N1,N2,P1)
%khai báo các tọa độ phân bố của trạm phát sóng BS
Xmax=P1(1); Xmin=P1(2); Ymax=P1(3); Ymin=P1(4); Zmax=P1(5); Zmin=P1(6);

N=N1*N2;
d=0.5;
%chạy x,y,z ngẫu nhiên trong khoảng đã cho
x = Xmax-rand*(Xmax-Xmin);
y = Ymax-rand*(Ymax-Ymin);
z = Zmax-rand*(Zmax-Zmin);
%hệ số độ lợi phức alpha theo phân phối chuẩn và được chuẩn hóa (tín hiệu tăng cường hoặc suy hao theo yếu tố môi trường ngẫu nhiên)
alpha = (normrnd(1, .1) + 1i*normrnd(1, .1)) / sqrt(2);
a = zeros(N,1);
for n1=1:N1
    for n2=1:N2
        a((n1-1)*N2+n2)=alpha*exp(-1j*2*pi*sqrt((x-(n1-1-(N1-1)/2)*d)^2+(z-(n2-1-(N2-1)/2)*d)^2+y^2));
        %tìm vector định hướng kênh BS-RIS bằng cách phức hóa modul vector từ điểm lấy mẫu BS tới phần tử RIS và nhân với hệ số độ lợi phức alpha
    end
end
G = 1*a;%gán giá trị vector a vào mảng G
%mỗi giá trị của mảng là vector định hướng kênh với biên độ là module vector và pha là arctan(ảo/thực)
