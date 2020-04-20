function P = medfuncParallebeamForwardProjection(theta, N, N_d, E)
% parallel beam forward projection function
% 输入参数说明：
% theta：投影角度矢量 in degrees
% N： 图像大小
% N_d：探测器通道个数
% E：仿真头部模型参数
% 输出说明：
% P：投影数据矩阵（N_d * theta_num)
shep = E;
theta_num = length(theta);
P = zeros(N_d,theta_num);
rho = shep(:,1);    %密度
ae = 0.5*N*shep(:,2);   %椭圆长轴
be = 0.5*N*shep(:,3);   %椭圆短轴
xe = 0.5*N*shep(:,4);   %x中心坐标
ye = 0.5*N*shep(:,5);   %y中心坐标
alpha = shep(:,6);  %对应椭圆旋转角度
alpha = alpha*pi/180;   %转换成弧度
theta = theta*pi/180;
TT = -(N_d-1)/2 : (N_d-1)/2;    %探测器坐标
for k1 = 1:theta_num
    p_theta = zeros(1,N_d);
    for k2 = 1:max(size(xe))
        rr = (ae(k2)*cos(theta(k1) - alpha(k2)))^2 + (be(k2)*sin(theta(k1) - alpha(k2)))^2;
        tt = (TT - xe(k2)*cos(theta(k1)) - ye(k2)*sin(theta(k1))).^2;
        temp = rr - tt;
        ind = temp>0;   %根号内值须为非负，这条语句意思时建立一个逻辑索引
        p_theta(ind) = p_theta(ind) + rho(k2)*(ae(k2)*be(k2)*sqrt(temp(ind)))./rr;
    end
    P(:,k1) = p_theta.';
end

