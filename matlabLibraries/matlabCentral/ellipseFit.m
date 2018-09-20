function report=ellipseFit(XY)
%ELLIPSEFIT - form 2D ellipse fit to given x,y data
%
%  report=ellipsefit(XY)
%
%in:
%
% XY: Input matrix of 2D coordinates to be fit. Each column XY(:,i) is [xi;yi]
%
%out: Finds the ellipse fitting the input data parametrized both as
%     A*x^2+B*x*y C*y^2+D*x+E*y=1 and [x-x0,y-y0]*Q*[x-x0;y-y0]=1
%
% report: a structure output with the following fields
%    
%    report.Q: the matrix Q
%    report.d: the vector [x0,y0]
%    report.ABCDE: the vector [A,B,C,D,E]
%    report.AxesDiams: The minor and major ellipse diameters
%    report.theta: The counter-clockwise rotation of the ellipse.
%
%NOTE: The code will give errors if the data fit traces out a non-elliptic or
%      degenerate conic section.
%
%See also ellipsoidfit
X=XY(1,:).';
Y=XY(2,:).';
M= [X.^2, X.*Y, Y.^2, X, Y, -ones(size(X,1),1)];
[U,S,V]=svd(M,0);
ABCDEF=V(:,end);
if size(ABCDEF,2)>1
      error 'Data cannot be fit with unique ellipse'
  else
      ABCDEF=num2cell(ABCDEF);
end
[A,B,C,D,E,F]=deal(ABCDEF{:});
Q=[A, B/2;B/2 C];
x0=-Q\[D;E]/2;
dd=F+x0'*Q*x0;
Q=Q/dd;
[R,eigen]=eig(Q);
eigen=eigen([1,4]);
if ~all(eigen>=0), %error 'Fit produced a non-elliptic conic section'; 
    report = 0;
else
    idx=eigen>0;
    eigen(idx)=1./eigen(idx);
    AxesDiams = 2*sqrt(eigen);
    theta=atand(tand(-atan2(R(1),R(2))*180/pi));
    report.Q=Q;
    report.d=x0(:).';
    report.ABCDE=[A, B, C, D, E]/F;
    report.AxesDiams=sort(AxesDiams(:)).';
    report.theta=theta;
end