function [ value ] = inter_sqaure( square_data,x,y )
%INTER_SQAURE �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
value = 0;
value = value+square_data(1,1)*(1-x)*(1-y);
value = value+square_data(1,2)*(1-x)*y;
value = value+square_data(2,1)*x*(1-y);
value = value+square_data(2,2)*x*y;

end

