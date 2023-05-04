%%% Descripción: Esta función se utiliza para realizar el analisis
%%% pseudo-esatico de un muro de retención  voladizo sección Cantilever  por el metodo de Simplificado.

%%% Entrdas: Base1: es la base 1 del muro de retención.  Base2: es la base 2 del muro de retención. Co:es la corona del
%%% muro. Ha:es la altura total del muro. D:es el ángulo delta, ángulo de fricción del muro.
%%% B: es el ángulo beta, ángulo de inclinación del muro. alpha: es el ángulo de inclinacíon del terreno.
%%% phib: es el ánsulo de fricción de la base del muro. concreto: es el
%%% peso específico del concreto. espesor: es la altura del estrato de suelo. 
%%% phi: es el ángulo de fricción de suelo de relleno. PesoE: es el peso específico del suelo de relleno.
%%% Cohesión: es la cohesión del suelo en la base del muro. n: es labcantidad de estratos.
%%% k: es el valor para identificar el estado del suelo. espe: es el espesor del suelo posterior al muro.
%%% ph: es el ángulo de fricción del suelo posterior.  R: es el valor para identificar la presencia de suelo posterior.
%%% TS: es el valor para identificar la presencia de agua en el suelo. R1 Y R2: son los valores para identificar el sistema de unidades.
%%% kh: es el coeficiente de aceleración horizontal del sismo. kv: es el coeficiente de aceleración vertical del sismo.
%%% So: es el valor para identificar la presencia de sobrecarga en el terreno.

%%% Salidas: FSD: es el factor de seguridad al desplazamiento sin tener encuenta el sismo.
%%% FSVR: es el factor de seguridad al volcamiento sumando a las fuerzas resistente el peso vertical sin tener encuenta el sismo.
%%% FSVA: es el factor de seguridad al volcamiento restando a las fuerzas actuantes el peso vertical sin tener encuenta el sismo.
%%% FSDs: es el factor de seguridad al desplazamiento teniendo encuenta el sismo.
%%% FSVRs: es el factor de seguridad al volcamiento sumando a las fuerzas resistente el peso vertical teniendo encuenta el sismo.
%%% FSVAs: es el factor de seguridad al volcamiento restando a las fuerzas actuantes el peso vertical teniendo encuenta el sismo.
%%% dEae: es el incremento sismico.
function [FSD,FSVR,FSVA,FSDs,FSVRs,FSVAs,dEae]= MtrizapseudoSIM(Base1,Base2,Base3,Ha,ht,tal,D,B,alpha,phib,concreto,espesor,phi,PesoE,Cohesion,n,k,espe,pe,ph,R,TS,R1,R2,kh,kv,as,So)
    if ismcc || isdeployed
        % Make sure DOM is compilable
        makeDOMCompilable();
     end
    import mlreportgen.dom.*
    import mlreportgen.report.*
    
    h=espesor;
    g=PesoE;
    c=Cohesion;   
    phii=phi;
    

    assignin('base','phib',phib);
    %Calculo de pie
    pie=Base3-Base2-Base1-tal;
    assignin('base','pie',pie);
    %Si hay suelo posterior
    assignin('base','R',R);
    if R==1
        a=0;
        nume=sind(B-ph)^2;
        p11=(sind(B)^2)*sind(B+D);
        p22=(1-sqrt(sind(ph+D)*sind(ph+a)/(sind(B+D)*sind(a+B))))^2;
        deno=p11*p22;
        kp=nume/deno;
        assignin('base','kp',kp);
        %Esfuerzo vertical
        ssv=espe*pe;
        assignin('base','espe',espe);
        assignin('base','pe',pe);
        %Esfuerzo horizontal
        ssh=ssv*kp;
        assignin('base','ssh',ssh);
        %Fuerza
        Fp=ssh*espe/2;
        Fps=Fp*sind(D);
        Fpc=Fp*cosd(D);
        %Brazo
        brap=espe/3;
        bras=Base2*2/3;
        brac=Ha/3;
        assignin('base','Fp',Fp);
        assignin('base','Fpc',Fpc);
        assignin('base','Fps',Fps);
        assignin('base','brap',brap);
        assignin('base','bras',bras);
        assignin('base','brac',brac);
        %Momento
        Mp=Fp*brap;
        assignin('base','Mp',Mp);
        
        Mps=Fps*bras;
        Mpc=Fpc*brac;

    end
    
    %Ka coulomb 
    %hallar ka o ko 
    if k==2    %cuando es activo
        num=sind(B+phii)^2;
        p1=(sind(B)^2)*sind(B-D);
        p2=(1+sqrt(sind(phii+D)*sind(phii-alpha)/(sind(B-D)*sind(alpha+B))))^2;
        den=p1*p2;
        ka=num/den;
        
        
    end
    if k==3 %cuando es pasivo
        nume=sind(B-phii)^2;
        p11=(sind(B)^2)*sind(B+D);
        p22=(1-sqrt(sind(phii+D)*sind(phii+a)/(sind(B+D)*sind(a+B))))^2;
        deno=p11*p22;
        ka=nume/deno;
        
    end

   
    
    %ACTUANTE TS TIPO DE SUELO SECO O CON AGUA
    if n==1
        if TS==3
            if R1==1
                Pactuante=1/2*(g-9.81)*h^2*ka;
            end
            if R2==1
                Pactuante=1/2*(g-62.4)*h^2*ka;
            end
            pvertical=Pactuante*sind(D);
            phor=Pactuante*cosd(D);
            if R1==1
                Pesoagua=1/2*9.81*h^2;
            end
            if R2==1
                Pesoagua=1/2*62.4*h^2;
            end
            
            %Brazo
            Ba=h/3;
            Bagua=h/3;
            Bvertical=Base3;
            assignin('base','Bvertical',Bvertical);
            %Momento
            Ma(1)=phor*Ba;
            Ma(2)=Pesoagua*Bagua;
            Mact=sum(Ma);
            Mvertical=Bvertical*pvertical;
            
            
            
        end
        if TS==2
            Pactuante=1/2*(g)*h^2*ka;
            pvertical=Pactuante*sind(D);
            phor=Pactuante*cosd(D);
            Pesoagua=0;
            
            %Brazo
            Ba=h/3;
            Bagua=h/3;
            Bvertical=Base3;
            %Momento
            Ma(1)=phor*Ba;
            Ma(2)=Pesoagua*Bagua;
            Mact=sum(Ma);
            Mvertical=Bvertical*pvertical;
        end
                
    end
    assignin('base','Bvertical',Bvertical);
    assignin('base','Pact',Pactuante);
    assignin('base','PVER',pvertical);
    assignin('base','PHOR',phor);
    assignin('base','Ma',Ma);
    assignin('base','Mact',Mact);
    assignin('base','Mvertical',Mvertical);

     %%%%%%%%%%%%%%%%%%%%%%SI HAY SOBRECARGA
    if So>0

        Soo=So*ka*(sind(B)^2/sind(B+alpha));
        Pcarga=Soo*h;
        Bcarga=h/2;
        Pecargav=Pcarga*sind(D);
        Pcargah=Pcarga*cosd(D);
        zt=(Pactuante*Ba+Pcarga*Bcarga)/(Pactuante+Pcarga);
        assignin('base','zt',zt);
    else
        Pcarga=0;
        zt=Ba;
    end

    Pactuante=Pactuante+Pcarga;
    pvertical=Pactuante*sind(D);
    phor=Pactuante*cosd(D);
    
    %Momento
    Ma(1)=phor*zt;
    Ma(2)=Pesoagua*Bagua;
    Mact=sum(Ma);
    Mvertical=Bvertical*pvertical;
    
    assignin('base','PCARG',Pcarga);
    assignin('base','pactuantec',Pactuante);
    assignin('base','phorc',phor);
    assignin('base','Mactc', Mact);
    assignin('base','Mverticalc', Mvertical);     
    assignin('base','pverc', pvertical);

    %RESISTENTE
    if alpha>0
        inclinada=tand(alpha)*tal;
        voli=inclinada*0.5*tal;        
        pinc=voli*g;
        assignin('base','Pincli',pinc);
        bi=tal*2/3+Base1+Base2+pie;
        momi=bi*pinc;
        assignin('base','momi',momi);
    else
        pinc=0;
        momi=0;
        assignin('base','momi',momi);
    end
     % COMO TIENE ZAPATA EL SUELO ACTUA TAMBIEN COMO RESISTENTE

    %Volumen PARTE DE CONCRETO
    Vol(1)=Base1*Ha*1;
    Vol(2)=Base2*Ha*1/2;
    Vol(3)=Base3*ht;
    %suelo
    Vol(4)=tal*Ha;
    %Peso
    Pes(1)=Vol(1)*concreto;%Kn
    Pes(2)=Vol(2)*concreto;%Kn
    Pes(3)=Vol(3)*concreto;%Kn
    %suelo
    Pes(4)=Vol(4)*g;%Kn
    assignin('base','PES',Pes);

    %Brazos
    Br(1)=Base2+Base1/2+pie;
    Br(2)=Base2*2/3+pie; 
    Br(3)=Base3/2; 
    Br(4)=tal/2+Base1+Base2+pie; 
    assignin('base','Br',Br);
    %Momentos
    M(1)=Br(1)*Pes(1);
    M(2)=Br(2)*Pes(2);
    M(3)=Br(3)*Pes(3);
    M(4)=Br(4)*Pes(4);
    assignin('base','M',M);
    Mr=sum(M)+momi;
    assignin('base','MR',Mr);


    %PARTE PSEUDOESTATICA
    
    teta=atand(kh/(1-kv));
    assignin('base','teta',teta);
    
   %METODO SIMPLIFICADO
%     kad=cosd(fitebe)^2/(cosd(teta)*(cosd(B)^2)*cosd(debeta)*dene);
    if kh>0 && kh<0.2
        dkad=(1.155-phi/75)*kh;
    end
    if kh>0.2 && kh<=0.35
        dkad=(phi/560)-0.08+((1.6-phi/45)*kh);
    end
    kad=ka+dkad;
    assignin('base','kad',kad);
    assignin('base','dkad',dkad);
    Ead=0.5*g*h^2*(1-kv)*kad;
        
        
  
    dEae=Ead-phor;
    assignin('base','Ead',Ead);
    assignin('base','dEae',dEae);
    %%%momento
    he=(((phor)*(h/3))+(dEae*0.6*h))/Ead;
    assignin('base','he',he);
    mEad=he*Ead*cosd(D);
    assignin('base','mEad',mEad);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Como se adiciono lo inclinado se tiene encuenta en el resistenre
    PesoRes=sum(Pes)+pinc;
   
    assignin('base','PesoRes',PesoRes);
    %FACTORES DE SEGURIDAD
    FSD=(((PesoRes+pvertical)*tand(phib))+(c*(Base1+Base2)))/(phor+ Pesoagua);
    FSVR=(Mr+Mvertical)/(Mact);
    FSVA=Mr/(Mact-Mvertical);
    assignin('base','FSVR',FSVR);
    
    %%%%%%%%%%%%%%%%%%%5
   

    %FACTORES DE SEGURIDAD sismo
    FSDs=(((PesoRes+pvertical)*tand(phib))+(c*(Base1+Base2)))/(phor+ Pesoagua+dEae);
    FSVRs=(Mr+Mvertical)/(Mact+mEad);
    FSVAs=Mr/(Mact-Mvertical+mEad);
    assignin('base','FSVRs',FSVRs);
    if R==1
        FSD=((PesoRes+pvertical)*tand(phib))/(phor+ Pesoagua-Fpc);
        FSVR=(Mr+Mvertical)/(Mact-Mps);
        FSVA=Mr/(Mact-Mps-Mvertical);
        
        FSDs=((PesoRes+pvertical)*tand(phib))/(phor+ Pesoagua-Fpc+dEae);
        FSVRs=(Mr+Mvertical)/(Mact-Mps+mEad);
        FSVAs=Mr/(Mact-Mps-Mvertical+mEad);
    end

      %%%%%%%%%%%%%%%GENERAR REPORTE

     reporte=questdlg('Desea generar reporte?', ...
	                        'Generar reporte', ...
	                        'SI',' NO  ','  CANCELAR  ','');
                        % Handle response
                        switch reporte
                            case 'SI' 
                                                             
           

                                % Add content to container (required)
                                % Types of content added here: title 
                                % page and table of contents reporters
                                doctype = "pdf";                  
                                rpt = Report("REPORTE MUROS DE RETENCIÓN",doctype);  
                                open(rpt);
                                
                    
                                layoutObj = rpt.Layout;
                                layoutObj.PageSize = PageSize("10in","11in","portrait");
                    
                                %Create Report container
                                rpt.Layout.Landscape = 1;
                                rpt.Layout.PageNumberFormat = "i";
                                titlepg = TitlePage;
                    
                                 %%%%%%%%%%%%
                                titlepg.Title = 'ANÁLISIS GEOTÉCNICO DE MUROS DE RETENCIÓN';
                                titlepg.Author = 'Juan Solarte';
                                tplayoutObj = titlepg.Layout;
                                tplayoutObj.PageSize = PageSize("10in","12in","portrait");
                                imDir = which("m.png");
                                if exist(imDir,"file") == 2
                                    ImObj = Image(imDir);
                                    ImObj.Height = "5in";
                                    ImObj.Width = "6in";
                                    titlepg.Image = ImObj;
                                end
                                add(rpt,titlepg);
                                add(rpt,TableOfContents);
                        
                    %%%%%%%%%%%%%%%5
                                % Add content to report sections (optional)
                                % Text and formal image added to chapter
                    
                    
                                chap = Chapter('Valores de entrada');
                                %TABLA MURO
                                s1=Section("CARACTERISTICAS MURO");
                                seccion1=["ALTURA MURO","BASE 1","BASE 2 ","BASE 3 " ,"ALTURA DE CUÑA (H.T)","PESO ESPECIFICO (concreto)","ANGULO DE INCLINACIÓN DEL MURO (β)"]';
                                volumen1=[Ha,Base1,Base2,Base3,ht,concreto,B]';
                               %Base1,Base2,Base3,Ha,ht,tal,D,B,alpha,phib,concreto,espesor,phi,PesoE,Cohesion,n,k,espe,pe,ph,R,TS,R1,R2,kh,kv,So
                                
                                Tabla1= FormalTable(["PROPIEDAD","VALOR"],[seccion1,volumen1]);
                                Tabla1.Header.Style{end+1} = BackgroundColor("silver");
                                Tabla1.Width = "500pt";
                                Tabla1.Border="solid";
                                Tabla1.ColSep="solid";
                                Tabla1.RowSep="solid";
                                Tabla1.HAlign="left";
                                Tabla1.TableEntriesHAlign="center";
                                
                                %TABLA SUELO RELLENO
                                 ss1=Section("CARACTERISTICAS SUELO DE RELLENO")
                                if TS==3
                                    seccions1=["PESO ESPECIFICO (suelo saturado)","ANGULO DE FRICCION (ϕ1)","ANGULO DE FRICCION ENTRE EL SUELO Y EL MURO (δ)","INCLINACION DEL RELLENO (α)","COHESIÓN","SOBRE-CARGA","COEFICIENTE HORIZONTAL DE ACERLERACIÓN(Kh)","COEFICIENTE VERTICAL DE ACERLERACIÓN(Kv)"]';
                                end
                                if TS==2
                                    seccions1=["PESO ESPECIFICO (suelo seco)","ANGULO DE FRICCION (ϕ1)","ANGULO DE FRICCION ENTRE EL SUELO Y EL MURO (δ)","INCLINACION DEL RELLENO (α)","COHESIÓN","SOBRE-CARGA","COEFICIENTE HORIZONTAL DE ACERLERACIÓN(Kh)","COEFICIENTE VERTICAL DE ACERLERACIÓN(Kv)"]';
                                end
                               
                                
                                volumens1=[PesoE,phi,D,alpha,0,So,kh,kv]';
                               
                                  %Base1,Base2,Base3,Ha,ht,tal,D,B,alpha,phib,concreto,espesor,phi,PesoE,Cohesion,n,k,espe,pe,ph,R,TS,R1,R2,kh,kv,So
                                
                                Tablas1= FormalTable(["PROPIEDAD","VALOR"],[seccions1,volumens1]);
                                Tablas1.Header.Style{end+1} = BackgroundColor("silver");
                                Tablas1.Width = "500pt";
                                Tablas1.Border="solid";
                                Tablas1.ColSep="solid";
                                Tablas1.RowSep="solid";
                                Tablas1.HAlign="left";
                                Tablas1.TableEntriesHAlign="center";
                    
                                %TABLA SUELO DE CIMENTACION
                                ssr1=Section("CARACTERISTICAS SUELO DE CIMENTACIÓN");
                                seccionsr1=["ANGULO DE FRICCION Base (ϕ2)","PESO ESPECIFICO (cimentación)","COHESIÓN"]';
                                volumensr1=[phib,PesoE,Cohesion]';
                               %phib,concreto,espesor,phi,PesoE,Cohesion,n,k,espe,pe,ph,R,TS,R1,R2,kh,kv,So
                                
                                
                                Tablasr1= FormalTable(["PROPIEDAD","VALOR"],[seccionsr1,volumensr1]);
                                Tablasr1.Header.Style{end+1} = BackgroundColor("silver");
                                Tablasr1.Width = "500pt";
                                Tablasr1.Border="solid";
                                Tablasr1.ColSep="solid";
                                Tablasr1.RowSep="solid";
                                Tablasr1.HAlign="left";
                                Tablasr1.TableEntriesHAlign="center";
                                
                                add(s1,Tabla1)
                                add(ss1,Tablas1)
                                add(ssr1,Tablasr1)
                                add(chap,s1)
                                add(chap,ss1)
                                add(chap,ssr1)
                                add(rpt,chap)
                                
                                 %%%%%%%Capitulo2
                                chap2 = Chapter('Factores Resistentes');
                                s=Section("Partes resistentes");
                                image2 = mlreportgen.report.FormalImage();
                                image2.Image = which('reportetriza.png');
                                image2.Caption = 'Partes Resistente(Muro Canteliver)';
                                image2.Height = '3in';
                                image2.Width = '4in';                             
                                                             
                                add(s,image2);
                               
%                                 plot1 = Image("reportetriza.png");
%                                 plot1.Width = "4in";
%                                 plot1.Height = "3in";
%                                 plot1.Style = {OuterMargin("150pt", "0pt","0pt","5pt")};
%                                 append(s,plot1);
                                
                                if alpha>0

                                    seccion=["A","B","C","D1","D2"," "]';
                                    volumen=[Vol(1),Vol(2),Vol(3),Vol(4),voli,Vol(1)+Vol(2)+Vol(3)+Vol(4)+voli]';
                                    peso=[Pes(1),Pes(2),Pes(3),Pes(4),pinc,PesoRes]';
                                    brazo=[round(Br(1),3),round(Br(2),3),round(Br(3),3),round(Br(4),3),round(bi,3)," "]';
                                    momento=[M(1),M(2),M(3),M(4),momi,Mr]';

                                    volumen=round(volumen,3);

                                    peso=round(peso,3);
                                    momento=round(momento,3);



                                    Tabla= FormalTable(["SECCIÓN","VOLUMEN","PESO","BRAZO","MOMENTO"],[seccion,volumen,peso,brazo,momento]);
                                    Tabla.Style = {OuterMargin("0pt", "0pt","50pt","5pt")};
                                    Tabla.Header.Style{end+1} = BackgroundColor("silver");
                                    Tabla.Width = "500pt";
                                    Tabla.Border="solid";
                                    Tabla.ColSep="solid";
                                    Tabla.RowSep="solid";
                                    Tabla.HAlign="center";
                                    Tabla.TableEntriesHAlign="center";

                                end

                                if alpha==0
                            
                                    seccion=["A","B","C","D"," "]';
                                    volumen=[Vol(1),Vol(2),Vol(3),Vol(4),Vol(1)+Vol(2)+Vol(3)+Vol(4)]';
                                    peso=[Pes(1),Pes(2),Pes(3),Pes(4),PesoRes]';
                                    brazo=[Br(1),Br(2),Br(3),Br(4)," "]';
                                    momento=[M(1),M(2),M(3),M(4),Mr]';
                                    Tabla= FormalTable(["SECCIÓN","VOLUMEN","PESO","BRAZO","MOMENTO"],[seccion,volumen,peso,brazo,momento]);
                                    Tabla.Style = {OuterMargin("0pt", "0pt","50pt","5pt")};
                                    Tabla.Header.Style{end+1} = BackgroundColor("silver");
                                    Tabla.Width = "500pt";
                                    Tabla.Border="solid";
                                    Tabla.ColSep="solid";
                                    Tabla.RowSep="solid";
                                    Tabla.HAlign="center";
                                    Tabla.TableEntriesHAlign="center";
                                end
                                
                                
                                add(s,Tabla);
                                add(chap2,s);
                                add(rpt,chap2);
                                
                                %%%%%%%Capitulo3
                                chap3 = Chapter('Factores Actuantes');
                                s3=Section("Condición Activa Coulomb");


%                                 plot3 = Image("kacou.png");
%                                 plot3.Style = {OuterMargin("100pt", "0pt","0pt","5pt")};
%                                 append(s3,plot3);
                    
                                condicion=["Ka"]';
                                valor=[ka]';
                               
                                
                                Tabla3= FormalTable(["CONDICIÓN", " "],[condicion,valor]);
                                Tabla3.Header.Style{end+1} = BackgroundColor("silver");
                                Tabla3.Style = {OuterMargin("0pt", "0pt","30pt","30pt")};
                                Tabla3.Width = "200pt";
                                Tabla3.Border="solid";
                                Tabla3.ColSep="solid";
                                Tabla3.RowSep="solid";
                                Tabla3.HAlign="center";
                                Tabla3.TableEntriesHAlign="center";
                                Tabla3.Style= [Tabla3.Style {NumberFormat("%.2f")}];
                                add(s3,Tabla3);
                    
                                %%%%TABLA FUERZAS
                                s32=Section("Fuerzas descompuestas");
                                secc=["Pa","Agua"]';
                                fue=[Pactuante,Pesoagua]';
                                
                                
                                Tabla31= FormalTable(["SECCIÓN"," Fuerzas "],[secc,fue]);
                                Tabla31.Header.Style{end+1} = BackgroundColor("silver");
                                Tabla31.Style = {OuterMargin("0pt", "0pt","40pt","40pt")};
                                Tabla31.Width = "400pt";
                                Tabla31.Border="solid";
                                Tabla31.ColSep="solid";
                                Tabla31.RowSep="solid";
                                Tabla31.HAlign="center";
                                Tabla31.TableEntriesHAlign="center";
                                Tabla31.Style= [Tabla31.Style {NumberFormat("%.2f")}];
                                add(s3,Tabla31);
                                %%%%FUERZAS DESCOMPUESTAS
                                secc1=["Pav (Pa*sen(δ))","Pah (Pa*cos(δ))","Agua"]';
                                fue=[pvertical,phor,Pesoagua]';
                                bra=[Bvertical,Ba,Ba]';
                                mom=[Mvertical,Ma(1),Ma(2)]';
                                
                                fue=round(fue,3);
                                bra=round(bra,3);
                                mom=round(mom,3);
                                Tabla32= FormalTable(["Componentes"," Fuerzas ","Brazos","Momento "],[secc1,fue,bra,mom]);
                                Tabla32.Header.Style{end+1} = BackgroundColor("silver");
                                Tabla32.Style = {OuterMargin("0pt", "0pt","40pt","40pt")};
                                Tabla32.Width = "500pt";
                                Tabla32.Border="solid";
                                Tabla32.ColSep="solid";
                                Tabla32.RowSep="solid";
                                Tabla32.HAlign="center";
                                Tabla32.TableEntriesHAlign="center";
                                Tabla32.Style= [Tabla32.Style {NumberFormat("%.2f")}];
                                add(s32,Tabla32);
                               
                    
                                add(chap3,s3);
                                add(chap3,s32)
                                add(rpt,chap3);
                                
                                 %%%%%%%Capitulo4
                                 
                    
                                chap4 = Chapter('Factores por sismo');
                                s4=Section('Metodo SIMPLIFICADO');
                               


                                image41 = mlreportgen.report.FormalImage();
                                image41.Image = which('kh2.png');
                                image41.Caption = 'Ecucaion para hallar Kh ';                                
                                add(s4,image41);
                                 
                                 if kh>0 && kh<0.2
                                      image42 = mlreportgen.report.FormalImage();
                                      image42.Image = which('menor.png');
                                      image42.Caption = 'Ecucaion para el Incrementos sismico ';
                                      add(s4,image42);
                                           
                                end
                                if kh>0.2 && kh<=0.35

                                     image42 = mlreportgen.report.FormalImage();
                                     image42.Image = which('mayor.png');
                                     image42.Caption = 'Ecucaion para el Incrementos sismico ';
                                     add(s4,image42);
                                    
                                end

                                image43 = mlreportgen.report.FormalImage();
                                image43.Image = which('kadsim.png');
                                image43.Caption = 'Ecucaion para el coeficiente Kad ';
                                add(s4,image43);
                              
                                 
                             
                                image44 = mlreportgen.report.FormalImage();
                                image44.Image = which('pae.png');
                                image44.Caption = 'Ecucaion para el coeficiente de presión Pae cuando hay sismo ';
                                add(s4,image44);
                    

                                image45 = mlreportgen.report.FormalImage();
                                image45.Image = which('z.png');
                                image45.Caption = 'Ecucaion para la ubicacion del Incremento sismico ';
                                add(s4,image45);
%                     
%                                 plot44 = Image("mea.png");
%                                 plot44.Style = {OuterMargin("0pt", "0pt","0pt","5pt")};
%                                 append(s4,plot44);
                                image46 = mlreportgen.report.FormalImage();
                                image46.Image = which('mea.png');
                                image46.Caption = 'Ecucaion para el momento generado por el Incremento sismico ';
                                add(s4,image46);
%                     
%                                 plot45 = Image("mea.png");
%                                 plot45.Style = {OuterMargin("0pt", "0pt","0pt","5pt")};
%                                 append(s4,plot45);
                                %%%%%TABLARESUMEN
                                te=[dkad]';
                                kaetabla=[kad]';
                                paetabla=[Ead]';
                                dpaetabla=[dEae]';
                                ztabla=[he]';
                                mea=[mEad]';
                                Tabla4= FormalTable(["ΔKad","Kae o Kad","Pae","ΔPae","Z","Mea"],[te,kaetabla,paetabla,dpaetabla,ztabla,mea]);
                                Tabla4.Style = {OuterMargin("0pt", "0pt","50pt","5pt")};
                                Tabla4.Header.Style{end+1} = BackgroundColor("silver");
                                Tabla4.Width = "400pt";
                                Tabla4.Border="solid";
                                Tabla4.ColSep="solid";
                                Tabla4.RowSep="solid";
                                Tabla4.HAlign="center";
                                Tabla4.TableEntriesHAlign="center";
                                Tabla4.Style= [Tabla4.Style {NumberFormat("%.2f")}];
                                add(s4,Tabla4);
                    
                    
                                add(chap4,s4);
                                add(rpt,chap4);
                     
                    
                              
                                
                                 %%%%%%%Capitulo5
                    
                    
                                chap5 = Chapter('Factores de seguridad estaticos');
                                s5=Section("Factor de seguridad al desplazamiento");
                                image5 = mlreportgen.report.FormalImage();
                                image5.Image = which('fsdes.png');
                                image5.Caption = 'Ecucaion para el Factor de seguridad al deslizamiento ';
                                add(s5,image5);
%                                 plot5 = Image("fsdes.png");
%                                 plot5.Style = {OuterMargin("100pt", "0pt","0pt","5pt")};
%                                 append(s5,plot5);
                                s52=Section("Factor de seguridad al volcamiento");
                                image52 = mlreportgen.report.FormalImage();
                                image52.Image = which('fsvol.png');
                                image5.Caption = 'Ecucaion para el Factor de seguridad al volcamiento ';
                                add(s52,image52);

%                                 plot52 = Image("fsvol.png");
%                                 plot52.Style = {OuterMargin("100pt", "0pt","0pt","5pt")};
%                                 append(s52,plot52);
%                                 
                                s53=Section("Resultados");
                                %TABLARESULTADOS
                                fsde=[FSD]';
                                fsv1=[FSVR]';
                                fsv2=[FSVA]';
                                
                               
                                Tabla5= FormalTable(["F.S DES","F.S.1 VOL","F.S.2 VOL"],[fsde,fsv1,fsv2]);
                                Tabla5.Style = {OuterMargin("0pt", "0pt","50pt","5pt")};
                                Tabla5.Header.Style{end+1} = BackgroundColor("silver");
                                Tabla5.Width = "400pt";
                                Tabla5.Border="solid";
                                Tabla5.ColSep="solid";
                                Tabla5.RowSep="solid";
                                Tabla5.HAlign="center";
                                Tabla5.TableEntriesHAlign="center";
                                Tabla5.Style= [Tabla5.Style {NumberFormat("%.2f")}];
                                add(s53,Tabla5);
                    
                    
                                add(chap5,s5);
                                add(chap5,s52);
                                add(chap5,s53);
                                add(rpt,chap5);
                                
                                 %%%%%%%Capitulo6
                    
                                chap6 = Chapter('Factores de seguridad pseudo-estaticos');
                                s6=Section("Factor de seguridad al desplazamiento por sismo");
                                image6 = mlreportgen.report.FormalImage();
                                image6.Image = which('fsdespe.png');
                                image6.Caption = 'Ecucaion para el Factor de seguridad al deslizamiento por sismo ';
                                add(s6,image6);
%                                 plot6 = Image("fsdespe.png");
%                                 plot6.Style = {OuterMargin("100pt", "0pt","0pt","5pt")};
%                                 append(s6,plot6);
                                s62=Section("Factor de seguridad al volcamiento por sismo");
%                                 plot62 = Image("fsvolpe.png");
%                                 plot62.Style = {OuterMargin("100pt", "0pt","0pt","5pt")};
%                                 append(s62,plot62);
                                image62 = mlreportgen.report.FormalImage();
                                image62.Image = which('fsvolpe.png');
                                image62.Caption = 'Ecucaion para el Factor de seguridad al volcamiento por sismo ';
                                add(s62,image62);

                                
                                s63=Section("Resultados");
                    
                                %TABLARESULTADOS
                                fsdepe=[FSDs]';
                                fsv1pe=[FSVRs]';
                                fsv2pe=[FSVAs]';
                               
                                
                               
                                Tabla6= FormalTable(["F.S DES (sismo)","F.S.1 VOL (sismo)","F.S.2 VOL(sismo)"],[fsdepe,fsv1pe,fsv2pe]);
                                Tabla6.Style = {OuterMargin("0pt", "0pt","50pt","5pt")};
                                 
                                Tabla6.Header.Style{end+1} = BackgroundColor("silver");
                                Tabla6.Width = "400pt";
                                Tabla6.Border="solid";
                                Tabla6.ColSep="solid";
                                Tabla6.RowSep="solid";
                                Tabla6.HAlign="center";
                                Tabla6.TableEntriesHAlign="center";
                                Tabla6.Style= [Tabla6.Style {NumberFormat("%.2f")}];
                                add(s63,Tabla6);
                              
                                add(chap6,s6);
                                add(chap6,s62);
                                add(chap6,s63);
                                add(rpt,chap6);
                    
                               
                                % Close the report (required)
                                close(rpt);
                                % Display the report (optional)
                                rptview(rpt);
                                

                                
                            case 'NO'
                                
                                
                            
                        end


    
    
    

end