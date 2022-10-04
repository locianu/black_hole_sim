import math
import numpy as np


def Sigma(r, theta):
    return pow(r, 2) + pow(a, 2) * pow(math.cos(theta), 2)
def drSigma(r):
    return 2 * r
def dthetaSigma(theta):
    return -2 * pow(a, 2) * math.cos(theta) * math.sin(theta)
def Delta(r):
    return pow(r, 2) - R * r + pow(a, 2)
def drDelta(r):
    return 2 * r - R

def main():

    plot_type = 1
    resolution_height = 2 #resoluções da tela, em pixels
    resolution_width = 2
    
    tStart = tic
    
    G = 1
    M = 1 #massa do buraco negro
    a = 0.6 #momento angular
    R = 2 * G * M
    radius_celestial_sphere = 80
    aDiskMin = 2 * R
    aDiskMax = 5 * R

#plot black sphere for the blackhole if plot_type = 1 is chosen and hold on the plot

    if plot_type == 1:
        [x, y, z] = sphere
        colormap([0,0,0])
        surf(R * x, R * y, R * z)
        axis equal
        hold on #isso aqui provavelmente não é feito assim em python
    
    #dimens�es da janela observacional
    window_height = 0.00001
    window_width = (resolution_width / resolution_height) * window_height
    distance_from_window = -1.4e-4

    coords_no_aDisk = np.zeros(resolution_height,resolution_width,3)
    coords_aDisk = np.zeros(resolution_height,resolution_width,3)
    
    stepsize = 0.1 #passo do rk4
    
    #hbar=parfor_progressbar(resolution_width,'Please wait ...') 
    #create the progress bar(resolution_width)
    
    fig = 1
    for j in range(0,resolution_width):
        #   hbar.iterate(1)
        #   update progress by one iteration
        for i in range(0,resolution_height):
            h = (window_height / 2) - (i - 1) * window_height / (resolution_height - 1)
            w = - (window_width / 2) + (j - 1) * window_width / (resolution_width - 1)
            r = 70
            theta = (math.pi / 2) - math.pi / 46 #offset the blackhole to see with of aDisk
            phi = 0
            t_dot = 1
            
            phi_dot = (w / math.sin(theta)) / math.sqrt((pow(a, 2) + pow(r, 2) * (pow(distance_from_window,2) + pow(w, 2) + pow(h, 2))))
            
            p_r = 2 * Sigma(r, theta) * (h * (pow(a, 2) + pow(r, 2)) * math.cos(theta) + r * sqrt(a^2+r^2) * math.sin(theta) * distance_from_window) / (math.sqrt(pow(distance_from_window, 2)+pow(h, 2)+pow(w, 2)) * (pow(a, 2) + 2 * pow(r, 2) + pow(a, 2) * math.cos(2 * theta)) * Delta(r))
            
            p_theta = 2 * Sigma(r,theta) * (-h * r * math.sin(theta) + math.sqrt(pow(a,2) + pow(r,2)) * math.cos(theta) * distance_from_window) / (math.sqrt(pow(distance_from_window, 2) + pow(h, 2) + pow(w, 2) * (pow(a, 2) + 2 * pow(r, 2) + pow(a, 2) * math.cos(2*theta))))
            E = (1 - R / r) * t_dot + (R * a * phi_dot) / r
            L = -( R * a) / r * t_dot + (pow(r, 2) + pow(a, 2) + (R * pow(a, 2)) / r) * phi_dot
            
            #geodesicas
            def f1(r,theta,p_r,p_theta):
                return (p_r * Delta(r)) / Sigma(r,theta)
            def f2(r,theta,p_r,p_theta):
                return (p_theta)/Sigma(r,theta)
            def f3(r,theta,p_r,p_theta):
                return -(1 / (2 * pow(Delta(r), 2) * pow(Sigma(r,theta), 2)) * (Sigma(r,theta) * (-E * Delta(r) * (a * R * (-2 * L + a * E * pow(math.sin(theta), 2)) + 2 * r * E * Sigma(r,theta)) + (a * (a * pow(L, 2) - 2 * L * r * R * E + a * r * R * pow(E, 2) * pow(math.sin(theta), 2)) + pow(p_r, 2) * pow(Delta(r), 2) + (pow(a, 2) + pow(r, 2)) * pow(E, 2) * Sigma(r,theta)) * drDelta(r)) + Delta(r) * (a * (L * (a * L - 2 * r * R * E) + a * r * R * pow(E, 2) * pow(math.sin(theta), 2)) - Delta(r) * (pow(p_theta, 2) + pow(L, 2) * (1 / pow(math.sin(theta), 2)) + pow(p_r, 2) * Delta(r))) * drSigma(r)))
            def f4(r,theta,p_r,p_theta):
                return -(1 / (2 * Delta(r) * pow(Sigma(r,theta), 2))) * (-2 * math.sin(theta) * (pow(a, 2) * r * R * pow(E, 2) * math.cos(theta) + pow(L, 2) * (1 / math.tan(theta)) * pow(1 / math.sin(theta), 3) * Delta(r)) * Sigma(r,theta) + (a * ( L * (a * L - 2 * r * R * E) + a * r * R * pow(E,2) * pow(math.sin(theta),2)) - Delta(r) * (pow(p_theta, 2) + pow(L, 2) * pow(1 / math.sin(theta), 2) + pow(p_r, 2) * Delta(r))) * dthetaSigma(theta))
            def f5(r,theta,p_r,p_theta):
                return (a * (-a * L + r * R * E) + L * pow(1 / math.sin(theta), 2) * Delta(r)) / (Delta(r) * Sigma(r,theta))
            
            x_0 = [r theta p_r p_theta phi]
            curve = x_0 #curve é uma função do matlab ou é só pra facilitar a chamada depois?
            
            k = 1
            Nk = 20000
            
            while ((R<r)&&(r<radius_celestial_sphere)&&(k<Nk)):
                #clean coordinates values
                curve(k, 2) = curve(k, 2) % (2 * math.pi)
                curve(k, 5) = curve(k, 5) % (2 * math.pi)
                if curve(k, 2) > math.pi:
                    curve(k,2) = 2*pi-curve(k,2)
                    curve(k,5) = (math.pi + curve(k,5)) % (2 * math.pi)
                theta = curve(k, 2)
                
                #runge-kutta
                step = min([stepsize * Delta(r), stepsize])
                k1 = step * f1(r,theta,p_r,p_theta)
                m1 = step * f2(r,theta,p_r,p_theta)
                n1 = step * f3(r,theta,p_r,p_theta) 
                s1 = step * f4(r,theta,p_r,p_theta)
                v1 = step * f5(r,theta,p_r,p_theta)
                
                k2 = step * f1(r+k1/2, theta+m1/2, p_r+n1/2, p_theta+s1/2)
                m2 = step * f2(r+k1/2, theta+m1/2, p_r+n1/2, p_theta+s1/2)
                n2 = step * f3(r+k1/2, theta+m1/2, p_r+n1/2, p_theta+s1/2)
                s2 = step * f4(r+k1/2, theta+m1/2, p_r+n1/2, p_theta+s1/2)
                v2 = step * f5(r+k1/2, theta+m1/2, p_r+n1/2, p_theta+s1/2)
                
                k3 = step * f1(r+k2/2, theta+m2/2, p_r+n2/2, p_theta+s2/2)
                m3 = step * f2(r+k2/2, theta+m2/2, p_r+n2/2, p_theta+s2/2)
                n3 = step * f3(r+k2/2, theta+m2/2, p_r+n2/2, p_theta+s2/2)
                s3 = step * f4(r+k2/2, theta+m2/2, p_r+n2/2, p_theta+s2/2)
                v3 = step * f5(r+k2/2, theta+m2/2, p_r+n2/2, p_theta+s2/2)
                
                k4 = step * f1(r+k3, theta+m3, p_r+n3, p_theta+s3)
                m4 = step * f2(r+k3, theta+m3, p_r+n3, p_theta+s3)
                n4 = step * f3(r+k3, theta+m3, p_r+n3, p_theta+s3)
                s4 = step * f4(r+k3, theta+m3, p_r+n3, p_theta+s3)
                v4 = step * f5(r+k3, theta+m3, p_r+n3, p_theta+s3)
                
                r = r + (k1 + (2*k2) + (2*k3) + k4)/6
                theta = theta + (m1 + (2 * m2) + (2 * m3) + m4)/6
                p_r = p_r + (n1 + (2 * n2) + (2 * n3) + n4)/6
                p_theta = p_theta + (s1 + (2 * s2) + (2 * s3) + s4)/6
                phi = phi + (v1 + (2 * v2) + (2 * v3) + v4)/6
                
                k = k+1
                disp(r)
                
                x = [r theta p_r p_theta phi]
                curve(k,:) = x #dúvida nisso
                
                #x=[r theta p_r p_theta];
                #disp(x);
                #disp(r);
                #disp(theta);
                #disp(["k=",k]);
            
            #Transform to euclidean coordinates and plot3
            [n,m]=size(curve)
            A=a*ones(n,1)
            
            PHIZ=ones(n,1)
    
            def Boyer2Cart(r,theta,phi):
                return [math.sqrt(r.^2+A.^2).*sin(theta).*cos(phi),math.sqrt(r.^2+A.^2).*sin(theta).*sin(phi),r.*cos(theta)]
                #porque os pontos nesse? é coisa do matlab?
            
            cart=Boyer2Cart(curve(:,1),curve(:,2),curve(:,5)) #porque os ':'?
    """%         disp([i,j]);
    %         pause;
          
    %         formatSpec = '%4.2f \t %4.2f \t %4.2f\n';
    %         filename = sprintf('plot%d.txt', fig);
    %         fileID = fopen(filename,'w');
    %         fprintf(fileID,formatSpec,cart(:,1),cart(:,2),cart(:,3));  % The format string is applied to each element of a
    %         fclose(fileID);
    %         fig = fig+1
    """        
           
            plot3(cart(:,1),cart(:,2),cart(:,3))
            
            """
    %         filename = sprintf('plot%d.png', fig);
    %         saveas(figi, filename);
    %         fig = fig+1;
             """       
"""    %close (hbar);
    % if plot_type == 1
    %     hold off
    % end
    
    %plot_r = curve(1,:)
    %plot_theta = curve(2,:)
    %polar(plot_theta,plot_r)
    %A = [plot_theta',plot_r']
    
    %disp(r);
    %disp(theta);"""
