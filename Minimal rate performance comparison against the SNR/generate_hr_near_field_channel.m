% định nghĩa hàm và khai báo các thông số đầu vào N1, N2, K và P
function [hK,x,y,z,alpha] = generate_hr_near_field_channel(N1,N2,K,P)

%trích xuất các giá trị tọa độ theo vector P
Xmax=P(1); Xmin=P(2); Ymax=P(3); Ymin=P(4); Zmax=P(5); Zmin=P(6);
%Xmax=0.2*P(1); Xmin=0.2*P(2); Ymax=0.2*P(3); Ymin=0.2*P(4); Zmax=0.2*P(5); Zmin=0.2*P(6);

N=N1*N2; %tính số phần tử mảng phản xạ
d=0.5;   %khoảng cách giữa các phần tử mảng phản xạ
hK=zeros(N,K); %khởi tạo ma trận hK với K cột và N hàng, mỗi hàng là 1 kênh truyền
aK=zeros(N,K); %khởi tạo ma trận aK với K cột và N hàng, mỗi hàng là 1 vector tín hiệu

% tạo vòng lặp cho từng người dùng
%%%% generate the channel h_{r,k} for each user k
for k=1:K  % biến chạy k chạy từ người dùng đầu tiên tới người dùng thứ K
    x = Xmax-rand*(Xmax-Xmin); % tạo giá trị x ngẫu nhiên chạy trong khoảng X
    y = Ymax-rand*(Ymax-Ymin); % tạo giá trị y ngẫu nhiên chạy trong khoảng Y
    z = Zmax-rand*(Zmax-Zmin); % tạo giá trị z ngẫu nhiên chạy trong khoảng Z

    % tạo phức alpha theo phân phối chuẩn
    alpha = (normrnd(1, .1) + 1i*normrnd(1, .1)) / sqrt(2);
    a = zeros(N,1); % khởi tạo biến a
    for n1=1:N1 % vòng lặp qua từng cột của ma trận kênh
        for n2=1:N2 % vòng lặp qua từng hàng của ma trận kênh
            %tính giá trị phần tử vector a với công thức kênh cận trường và hệ số alpha
            a((n1-1)*N2+n2)=alpha*exp(-1j*2*pi*sqrt((x-(n1-1-(N1-1)/2)*d)^2+(z-(n2-1-(N2-1)/2)*d)^2+y^2));
        end
    end
    hr = 1*a; % gán giá trị vector a vào kênh hr
    hK(:,k) = hr./sqrt(N); %gán giá trị vector kênh hr vào cột thứ K của ma trận kênh hK và chuẩn hóa theo căn bậc 2 của N
    aK(:,k) = a; % lưu vector tín hiệu a vào ma trận aK cho mỗi người dùng
end

