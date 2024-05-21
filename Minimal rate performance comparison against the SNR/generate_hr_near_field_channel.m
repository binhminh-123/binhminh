function [hK,x,y,z,alpha] = generate_hr_near_field_channel(N1,N2,K,P)
%khai báo các giá trị tọa độ phân bố của UE
Xmax=P(1); Xmin=P(2); Ymax=P(3); Ymin=P(4); Zmax=P(5); Zmin=P(6);

N=N1*N2;
d=0.5;
hK=zeros(N,K);
aK=zeros(N,K);

%%%% tạo kênh truyền cho từng người dùng
for k=1:K
    %các giá trị tọa độ x,y,z chạy ngẫu nhiên trong khoảng đã cho
    x = Xmax-rand*(Xmax-Xmin);
    y = Ymax-rand*(Ymax-Ymin);
    z = Zmax-rand*(Zmax-Zmin);

    alpha = (normrnd(1, .1) + 1i*normrnd(1, .1)) / sqrt(2);%hệ số độ lợi phức alpha theo phân phối chuẩn chuẩn hóa làm yếu tố môi trường
    a = zeros(N,1);
    for n1=1:N1
        for n2=1:N2
            a((n1-1)*N2+n2)=alpha*exp(-1j*2*pi*sqrt((x-(n1-1-(N1-1)/2)*d)^2+(z-(n2-1-(N2-1)/2)*d)^2+y^2));
            %tìm vector định hướng kênh RIS-UE bằng cách phức hóa modul vector từ phần tử RIS tới điểm lấy mẫu UE và nhân với hệ số độ lợi phức alpha
        end
    end
    hr = 1*a;%lưu giá trị vector a vào mảng kênh hr
    hK(:,k) = hr./sqrt(N);%lưu giá trị vector đã chuẩn hóa của mảng kênh hr vào cột thứ k mảng kênh hK
    aK(:,k) = a;%???
end

