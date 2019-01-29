#The formula for calculating tree likelihood is equation (4) in Gavryushkina et al (2014) multiplied by (m+M+k+K)!*q1(x1)/(lambda1*(1-p0hat(x1))^2)
lambda1<-1.5
lambda2<-0.5
mu1<-0.4
mu2<-0.3
psi1<-0.3
psi2<-0.8
rho<-0.9
A2<-sqrt((lambda2-mu2-psi2)^2+4*lambda2*psi2)
B2<-((1-2*(1-rho))*lambda2+mu2+psi2)/A2
p2from1<-(lambda2+mu2+psi2-A2*(exp(A2)*(1+B2)-(1-B2))/(exp(A2)*(1+B2)+(1-B2)))/(2*lambda2)
p2from0point8<-(lambda2+mu2+psi2-A2*(exp(0.8*A2)*(1+B2)-(1-B2))/(exp(0.8*A2)*(1+B2)+(1-B2)))/(2*lambda2)
q2from0point8<-4*exp(0.8*A2)/(exp(0.8*A2)*(1+B2)+(1-B2))^2
A1<-sqrt((lambda1-mu1-psi1)^2+4*lambda1*psi1)
B1<-((1-2*p2from1)*lambda1+mu1+psi1)/A1
q1from2point5<-4*exp(1.5*A1)/(exp(1.5*A1)*(1+B1)+(1-B1))^2
q2from1<-4*exp(A2)/(exp(A2)*(1+B2)+(1-B2))^2
A2hat<-sqrt((lambda2-mu2)^2)
B2hat<-((1-2*(1-rho))*lambda2+mu2)/A2hat
A1hat<-sqrt((lambda1-mu1)^2)
p2hatfrom1<-(lambda2+mu2-A2hat*(exp(A2hat)*(1+B2hat)-(1-B2hat))/(exp(A2hat)*(1+B2hat)+(1-B2hat)))/(2*lambda2)
B1hat<-((1-2*p2hatfrom1)*lambda1+mu1)/A1hat
p1hatfrom2point5<-(lambda1+mu1-A1hat*(exp(1.5*A1hat)*(1+B1hat)-(1-B1hat))/(exp(1.5*A1hat)*(1+B1hat)+(1-B1hat)))/(2*lambda1)
logFfromT<-log(2*psi1*psi2*rho*(q1from2point5)^2*p2from0point8) -2*log(1-p1hatfrom2point5)-log(q2from0point8)+2*log(q2from1)
print(logFfromT, digits=15)

              
