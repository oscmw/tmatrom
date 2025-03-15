classdef obstacleCustomShape < obstacle
    
    properties
        amplitude = 0.3; % Şeklin çıkıntılarının yüksekliği
        frequency = 4;   % Şeklin kenarlarının sayısı
    end
    
    methods
        function self = obstacleCustomShape()
            self = self@obstacle();
        end
        
        function [x, y, qx, qy, dqx, dqy, ddqx, ddqy, nqx, nqy, jac] = geom(self, t)
            % Çember koordinatları
            x = cos(t);
            y = sin(t);
            
            % Asimetrik şekli oluşturmak için radyal modifikasyon
            r = 1 + self.amplitude * sin(self.frequency * t) + 0.1 * cos(2 * t);
            qx = r .* cos(t);
            qy = r .* sin(t);
            
            % Türevler
            dr = self.amplitude * self.frequency * cos(self.frequency * t) - 0.2 * sin(2 * t);
            dqx = dr .* cos(t) - r .* sin(t);
            dqy = dr .* sin(t) + r .* cos(t);
            ddqx = -self.amplitude * self.frequency^2 * sin(self.frequency * t) .* cos(t) ...
                   - 2 * dr .* sin(t) - r .* cos(t);
            ddqy = -self.amplitude * self.frequency^2 * sin(self.frequency * t) .* sin(t) ...
                   + 2 * dr .* cos(t) - r .* sin(t);
            
            % Jacobian ve normal vektörler
            jac = sqrt(dqx.^2 + dqy.^2);
            nqx = dqy ./ jac;
            nqy = -dqx ./ jac;
        end
    end
end
