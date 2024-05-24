function [codebook,record] = generate_near_field_codebook(N1,N2,d,P1,P2,Delta)
%codebook được tạo từ giá trị của các vector định hướng, mỗi vector đc tạo từ 1 cặp điểm lấy mẫu trong 2 vùng phân bố
N=N1*N2;

Xmax1=P1(1); Xmin1=P1(2); Ymax1=P1(3); Ymin1=P1(4); Zmax1=P1(5); Zmin1=P1(6);%giới hạn vùng phân bố của trạm phát sóng
Xmax2=P2(1); Xmin2=P2(2); Ymax2=P2(3); Ymin2=P2(4); Zmax2=P2(5); Zmin2=P2(6);%giới hạn vùng phân bố của User
Xdelta1=Delta(1);Ydelta1=Delta(2);Zdelta1=Delta(3);Xdelta2=Delta(4);Ydelta2=Delta(5);Zdelta2=Delta(6);%khoảng cách giữa các điểm lấy mẫu trên lưới tọa độ

% Xgrid1 = linspace(Xmin1,Xmax1,Xnum1); Ygrid1 = linspace(Ymin1,Ymax1,Ynum1); Zgrid1 = linspace(Zmin1,Zmax1,Znum1);
% Xgrid2 = linspace(Xmin2,Xmax2,Xnum2); Ygrid2 = linspace(Ymin2,Ymax2,Ynum2); Zgrid2 = linspace(Zmin2,Zmax2,Znum2);
%tăng A -> Delta1 tăng -> khoảng nhảy tăng -> số mã giảm
%tọa độ grid sẽ chạy từ Xmin1 tới Xmax1 và chấm sau mỗi khoảng nhảy Xdelta1, Y và Z tương tự
Xgrid1=[Xmin1:Xdelta1:Xmax1]; Ygrid1=[Ymin1:Ydelta1:Ymax1]; Zgrid1=[Zmin1:Zdelta1:Zmax1];% tạo lưới tọa độ các điểm lấy mẫu cho kênh BS-RIS
Xgrid2=[Xmin2:Xdelta2:Xmax2]; Ygrid2=[Ymin2:Ydelta2:Ymax2]; Zgrid2=[Zmin2:Zdelta2:Zmax2];% tạo lưới tọa độ các điểm lấy mẫu cho kênh RIS-UE

Xnum1=length(Xgrid1); Ynum1=length(Ygrid1); Znum1=length(Zgrid1);%tính số điểm lưới của lưới kênh BS-RIS
Xnum2=length(Xgrid2); Ynum2=length(Ygrid2); Znum2=length(Zgrid2);%tính số điểm lưới của lưới kênh RIS-UE

record=zeros(Xnum1*Ynum1*Znum1*Xnum2*Ynum2*Znum2,6);%tạo mảng lưu giá trị 2 điểm trên 2 lưới kênh
codebook = zeros(Xnum1*Ynum1*Znum1*Xnum2*Ynum2*Znum2,N);%tạo ma trận lưu các vector mã codeword, phép nhân để xác định tổng số điểm mẫu cần thiết
i=1;
%vòng lặp quét qua từng điểm của lưới kênh BS-RIS
for x1=Xgrid1
    for y1=Ygrid1
        for z1=Zgrid1
            %vòng lặp quét qua từng điểm của lưới kênh RIS-UE
            for x2=Xgrid2
                for y2=Ygrid2
                    for z2=Zgrid2
                        a = zeros(1,N);
                        for n1=1:N1
                            for n2=1:N2
                                a((n1-1)*N2+n2) = exp(1j*2*pi*(sqrt((x1-(n1-1-(N1-1)/2)*d)^2+(z1-(n2-1-(N2-1)/2)*d)^2+y1^2)+sqrt((x2-(n1-1-(N1-1)/2)*d)^2+(z2-(n2-1-(N2-1)/2)*d)^2+y2^2)));
  %tính vector định hướng bằng cách tìm 2 vector từ 2 điểm lưới tới phần tử của mảng RIS, sau đó tính tổng module sẽ được độ dài kênh truyền, lấy phức hợp ta sẽ ra đc vector định hướng tương ứng chứa pha và biên độ tín hiệu phù hợp với kênh đó. x1-(n1-1-(N1-1)/2)*d trung tâm hóa tọa độ phần tử phản xạ, lý do chọn tọa độ theo công thức đó là để các 2 điểm được chọn luôn nằm đối xứng qua 0.                           
                            end
                        end
%                         a=a/sqrt(N);
                        codebook(i,:)=a;%lưu giá trị vector a vào hàng thứ i trong tổng số điểm mẫu
                        record(i,:)=[x1,y1,z1,x2,y2,z2];%lưu giá trị tọa độ các điểm mẫu qua các vòng lặp
                        i=i+1;
                    end
                end
            end
        end
    end
end

[codebook,index]=unique(codebook,'row');%loại các hàng có giá trị trùng trong quá trình lặp
record=record(index,:);%cập nhật lại mảng record sau khi loại các giá trị trùng lặp
end

