//JavaScript code for simulation of neutron Laue diffraction pattern at HRC

//2020/6/18,  introduce lattife constants and sample orientation 
//2020/6/5
var version = "0.2";

var TOFconst = 2.286;       // TOF at 1 m is 2.286/sqrt(E)

var X0 = 0;
var Y0 = 250;

var length1=200;
var radius=5;

var HD = 40;    // height of center of PSD from incident beam
var LD = 2800;  // length of PSD
var L20 = 4000; // distance from sample to PSD in horizontal plane

var u = new Array(3);
var v = new Array(3);

var ux = new Array(3);
var vx = new Array(3);

var a_star = new Array(3);
var b_star = new Array(3);
var c_star = new Array(3);

a_star[0]=1.0;
a_star[1]=0.0;
a_star[2]=0.0;

b_star[0]=0.0;
b_star[1]=1.0;
b_star[2]=0.0;

c_star[0]=0.0;
c_star[1]=0.0;
c_star[2]=1.0;

var Hmax;
var Kmax;
var Lmax;

var maxphih = 60.0;
var maxphiv = 60.0;
var scaleX=800;
var scaleY=500;

var maxphih = 60.0;
var maxphiv = 60.0;
var scaleX=800;
var scaleY=500;


function draw() {
    document.getElementById("verNum").innerHTML=version;
    document.getElementById("verNum2").innerHTML=version;


    //set_Lattice();
    set_Lattice();
    draw_DetMap();

}

function set_Lattice(){

    //input parameters: lattice constants and sample orientation)
    a = Number(document.getElementById('a').value);
    b = Number(document.getElementById('a').value);
    c = Number(document.getElementById('a').value);
    alpha = Number(document.getElementById('alpha').value)/180.0*Math.PI;   // in radian
    beta  = Number(document.getElementById('beta').value)/180.0*Math.PI;    // in radian
    gamma = Number(document.getElementById('gamma').value)/180.0*Math.PI;   // in radian
    u[0] = Number(document.getElementById('u1').value);
    u[1] = Number(document.getElementById('u2').value);
    u[2] = Number(document.getElementById('u3').value);
    v[0] = Number(document.getElementById('v1').value);
    v[1] = Number(document.getElementById('v2').value);
    v[2] = Number(document.getElementById('v3').value);

    // calculation
    let DD = (Math.cos(alpha)-Math.cos(gamma)*Math.cos(beta))/Math.sin(gamma);
    let PP = Math.sqrt(Math.sin(beta)-DD**2.0);

    ux[0] = 2.0*Math.PI*u[0]/a;
    ux[1] = 2.0*Math.PI*(-u[0]/a/Math.tan(gamma)+u[1]/b/Math.sin(gamma));
    ux[2] = 2.0*Math.PI*(u[0]/a*(DD/Math.tan(gamma)-Math.cos(beta))-u[1]/b*DD/Math.sin(gamma)+u[2]/c)/PP;
    vx[0] = 2.0*Math.PI*v[0]/a;
    vx[1] = 2.0*Math.PI*(-v[0]/a/Math.tan(gamma)+v[1]/b/Math.sin(gamma));
    vx[2] = 2.0*Math.PI*(v[0]/a*(DD/Math.tan(gamma)-Math.cos(beta))-v[1]/b*DD/Math.sin(gamma)+v[2]/c)/PP;

    let uy2uz2 = ux[1]**2.0+ux[2]**2.0;    
    let Uabs = Math.sqrt(ux[0]**2.0+uy2uz2);

    let Rvy =(-vx[0]*ux[1]+(vx[1]*(ux[0]*ux[1]**2.0+Ubas*ux[2]**2.0)+vx[2]*ux[1]*ux[2]*(ux[0]-Ubas))/uy2uz2)/Uabs;
    let Rvz =(-vx[0]*ux[2]+(vx[2]*(ux[0]*ux[2]**2.0+Ubas*ux[1]**2.0)+vx[1]*ux[2]*ux[1]*(ux[0]-Ubas))/uy2uz2)/Uabs;

    let cosphi=Rvy/Math.sqrt(Rvy**2.0+Rvz**2.0);
    let sinphi=Rvz/Math.sqrt(Rvy**2.0+Rvz**2.0);



}

function draw_DetMap(){

    var canvas = document.getElementById('CanvasDetMap');
    var context = canvas.getContext('2d');

    //refresh
    context.clearRect(0, 0, canvas.width, canvas.height);
    context.strokeStyle = "rgb(0, 0, 0)";
    context.lineWidth=1;

    Hmax = Number(document.getElementById('Hmax').value);
    Kmax = Number(document.getElementById('Kmax').value);
    Lmax = Number(document.getElementById('Lmax').value);

    // line
    context.strokeStyle = "rgb(255, 0, 0)";
    context.beginPath();
    context.moveTo(X0, Y0);
    context.lineTo(X0+length1, Y0);
    context.stroke();

    // circle
    context.strokeStyle = "rgb(0, 150, 0)";
    var delta=20;
    var limit=5;

    var temp2=0;
    var Ghkl=new Array(3);
    for (var H=-Hmax;H<=Hmax;H+=1){
        for (var K=-Kmax;K<=Kmax;K+=1){
            for (var L=-Lmax;L<=Lmax;L+=1){

                if((H==0)&&(K==0)&&(L==0)){

                }
                else{
                    for(let i=0;i<3;i+=1){
                        Ghkl[i]=H*a_star[i]+K*b_star[i]+L*c_star[i];
                    }
    
                    let G_len=0;        //calculate length of G
                    for(let i=0;i<3;i+=1){
                        G_len=G_len+Ghkl[i]**2.0;
                    }
                    G_len=Math.sqrt(G_len);
                 
                    let sinphiv=Ghkl[2]/G_len;
                    let phiv = Math.asin(sinphiv);      // in radian
                    let sinphih=Ghkl[1]/(G_len*Math.cos(phiv));
                    let phih=  Math.asin(sinphih);      // in radian
                    //let cos2th=Ghkl[0]/G_len;
                    //let twotheta= Math.acos(cos2th);   //in radidan

                    let PosX=scaleX*phih/Math.PI*180.0/maxphih+X0;
                    //let PosY=scaleY*phiv/Math.PI*180.0/maxphiv+Y0;
                    let PosY=scaleY*(HD+LD/2-L20*Math.tan(phiv))/LD

                    context.beginPath();
                    context.arc(PosX,PosY, radius, 0, 2 * Math.PI);
                    context.stroke();    
                }
            }
        }
    }

    //text for debug

    context.font = "italic 13px sans-serif";
    //context.fillText(Ghkl[2], X0, Y0+length1);
    //context.fillText(alpha, X0, Y0+length1);
    // console.log("Hallo")
    //context.font = "italic 10px sans-serif";
    //context.fillText(Ghkl[2], X0, Y0+length1);


}

