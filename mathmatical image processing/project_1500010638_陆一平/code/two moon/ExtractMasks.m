function [DFilters RFilters]=ExtractMasks(type)
switch type
    case 'Haar'
        DFilters{1} = @(x) (cos(x/2));
        DFilters{2} = @(x) (sin(x/2));
        RFilters{1} = @(x) (cos(x/2));
        RFilters{2} = @(x) (sin(x/2));
    case 'Linear'
        DFilters{1} = @(x) (cos(x/2).^2);
        DFilters{2} = @(x) (sin(x)/sqrt(2));
        DFilters{3} = @(x) (sin(x/2).^2);
        RFilters{1} = @(x) (cos(x/2).^2);
        RFilters{2} = @(x) (sin(x)/sqrt(2));
        RFilters{3} = @(x) (sin(x/2).^2);
    case 'Cubic'
        DFilters{1} = @(x) (cos(x/2).^3);
        DFilters{2} = @(x) (sqrt(3)*sin(x/2).*cos(x/2).^2);
        DFilters{3} = @(x) (sqrt(3)*sin(x/2).^2.*cos(x/2));
        DFilters{4} = @(x) (sin(x/2).^3);
        RFilters{1} = @(x) (cos(x/2).^3);
        RFilters{2} = @(x) (sqrt(3)*sin(x/2).*cos(x/2).^2);
        RFilters{3} = @(x) (sqrt(3)*sin(x/2).^2.*cos(x/2));
        RFilters{4} = @(x) (sin(x/2).^3);
%         DFilters{1} = @(x) (exp(1i*x/2).*cos(x/2).^3);
%         DFilters{2} = @(x) (1i*exp(1i*x/2).*sqrt(3).*sin(x/2).*cos(x/2).^2);
%         DFilters{3} = @(x) (exp(1i*x/2).*sqrt(3).*sin(x/2).^2.*cos(x/2));
%         DFilters{4} = @(x) (-1i*exp(1i*x/2).*sin(x/2).^3);
%         RFilters{1} = @(x) (exp(-1i*x/2).*cos(x/2).^3);
%         RFilters{2} = @(x) (-1i*exp(-1i*x/2).*sqrt(3).*sin(x/2).*cos(x/2).^2);
%         RFilters{3} = @(x) (exp(-1i*x/2).*sqrt(3).*sin(x/2).^2.*cos(x/2));
%         RFilters{4} = @(x) (1i*exp(-1i*x/2).*sin(x/2).^3);
    case 'Pseudo-Spline31'
        DFilters{1} = @(x) (cos(x/2).^6.*(1+3*sin(x/2).^2));
        DFilters{2} = @(x) (exp(1i*x).*sin(x/2).^6.*(1+3*cos(x/2).^2));
        DFilters{3} = @(x) (0.5*((4 - cos(x/2).^12.*(3*cos(x) - 5).^2 - sin(x/2).^12.*(3*cos(x) + 5).^2).^(1/2)./2)+ ...
            exp(1i*x).*0.5.*((4 - cos(x/2).^12.*(3*cos(x) - 5).^2 - sin(x/2).^12.*(3*cos(x) + 5).^2).^(1/2)./2));
        DFilters{4} = @(x) (exp(1i*x).*0.5.*((4 - cos(x/2).^12.*(3*cos(x) - 5).^2 - sin(x/2).^12.*(3*cos(x) + 5).^2).^(1/2)./2)- ...
            0.5.*((4 - cos(x/2).^12.*(3*cos(x) - 5).^2 - sin(x/2).^12.*(3*cos(x) + 5).^2).^(1/2)./2));
        RFilters{1} = @(x) (cos(x/2).^6.*(1+3*sin(x/2).^2));
        RFilters{2} = @(x) (exp(-1i*x).*sin(x/2).^6.*(1+3*cos(x/2).^2));
        RFilters{3} = @(x) (0.5*((4 - cos(x/2).^12.*(3*cos(x) - 5).^2 - sin(x/2).^12.*(3*cos(x) + 5).^2).^(1/2)./2)+ ...
            exp(-1i*x).*0.5.*((4 - cos(x/2).^12.*(3*cos(x) - 5).^2 - sin(x/2).^12.*(3*cos(x) + 5).^2).^(1/2)./2));
        RFilters{4} = @(x) (exp(-1i*x).*0.5.*((4 - cos(x/2).^12.*(3*cos(x) - 5).^2 - sin(x/2).^12.*(3*cos(x) + 5).^2).^(1/2)./2)- ...
            0.5.*((4 - cos(x/2).^12.*(3*cos(x) - 5).^2 - sin(x/2).^12.*(3*cos(x) + 5).^2).^(1/2)./2));
end