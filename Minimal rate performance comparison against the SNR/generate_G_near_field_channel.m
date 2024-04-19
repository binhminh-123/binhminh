%định nghĩa hàm và các thông số đầu vào N1, N2 và P
function [G,x,y,z,alpha] = generate_G_near_field_channel(N1,N2,P)

%trích xuất các giá trị tọa độ từ vector P
Xmax=P(1); Xmin=P(2); Ymax=P(3); Ymin=P(4); Zmax=P(5); Zmin=P(6);

N=N1*N2; %tính số phần tử trong mảng phản xạ
d=0.5; %khoảng cách giữa các phần tử mảng

x = Xmax-rand*(Xmax-Xmin); %tạo giá trị x ngẫu nhiên chạy trong khoảng X
y = Ymax-rand*(Ymax-Ymin); %tạo giá trị y ngẫu nhiên chạy trong khoảng Y
z = Zmax-rand*(Zmax-Zmin); %tạo giá trị z ngẫu nhiên chạy trong khoảng Z

%tạo phức alpha với phần thực và phần ảo theo phân phối chuẩn
alpha = (normrnd(1, .1) + 1i*normrnd(1, .1)) / sqrt(2);
a = zeros(N,1); %khởi tạo biến vector a

for n1=1:N1 %vòng lặp qua từng cột của ma trận kênh truyền
    for n2=1:N2 %vòng lặp qua từng dòng của ma trận kênh truyền
                %tính giá trị mỗi phần tử của ma trận kênh G với công thức kênh cận trường và hệ số alpha
        a((n1-1)*N2+n2)=alpha*exp(-1j*2*pi*sqrt((x-(n1-1-(N1-1)/2)*d)^2+(z-(n2-1-(N2-1)/2)*d)^2+y^2));
    end
end
G = 1*a; %gán giá trị ma trận kênh G bằng giá trị vector a

