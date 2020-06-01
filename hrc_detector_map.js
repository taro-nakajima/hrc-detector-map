//JavaScript code for simulation of neutron Laue diffraction pattern at HRC
//test
var version = "0.1";
var TOFconst = 2.286;       // TOF at 1 m is 2.286/sqrt(E)

var X0 = 0;
var Y0 = 250;

var length1=200;
var radius=5;

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

var Hmax=5;
var Kmax=5;
var Lmax=1;

var maxphih = 60.0;
var maxphiv = 60.0;
var scaleX=800;
var scaleY=500;


function draw() {
    document.getElementById("verNum").innerHTML=version;
    document.getElementById("verNum2").innerHTML=version;

    draw_DetMap();

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
                    let phiv = Math.asin(sinphiv);
                    let sinphih=Ghkl[1]/(G_len*Math.cos(phiv));
                    let phih=  Math.asin(sinphih);
                    //let cos2th=Ghkl[0]/G_len;
                    //let twotheta= Math.acos(cos2th);   //in radidan

                    let PosX=scaleX*phih/Math.PI*180.0/maxphih+X0;
                    let PosY=scaleY*phiv/Math.PI*180.0/maxphiv+Y0;
    
                    context.beginPath();
                    context.arc(PosX,PosY, radius, 0, 2 * Math.PI);
                    context.stroke();    
                }
            }
        }
    }

    //text for debug
    //context.font = "italic 10px sans-serif";
    //context.fillText(Ghkl[2], X0, Y0+length1);
    

}

