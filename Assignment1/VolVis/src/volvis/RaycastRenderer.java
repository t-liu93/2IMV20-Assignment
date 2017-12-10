/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import gui.RaycastRendererPanel;
import gui.TransferFunction2DEditor;
import gui.TransferFunctionEditor;
import java.awt.image.BufferedImage;
import util.TFChangeListener;
import util.VectorMath;
import volume.GradientVolume;
import volume.Volume;
import volume.VoxelGradient;

/**
 *
 * @author michel
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {

    private Volume volume = null;
    private GradientVolume gradients = null;
    private boolean shadeOn = false;
    private static final double Ia = 0.1;
    private static final double kdiff = 0.7;
    private static final double kspec = 0.2;
    private static final int a = 10;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;
    
    //Variables for listeners
    private String status;
    
    //Control variables for MIP
    private int MIPSlices = 50; //The number of slices will be checked in MIP mode apart from the center slice
    //More slices will give a more precise visualization, but it is slower.
    //Dafault 50
    private int MIPSliceStep = 1; //Each time when checking another slice, the step to be applied on the previous one
    //The small step may give a more precise visualization but it is very slow. 
    //Default 1
    
    //Constructors    
    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
        
        status = "slicer";
    }

    public void setVolume(Volume vol) {
        System.out.println("Assigning volume");
        volume = vol;

        System.out.println("Computing gradients");
        gradients = new GradientVolume(vol);

        // set up image for storing the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ()));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
        // create a standard TF where lowest intensity maps to black, the highest to white, and opacity increases
        // linearly from 0.0 to 1.0 over the intensity range
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());
        
        // uncomment this to initialize the TF with good starting values for the orange dataset 
        //tFunc.setTestFunc();
        
        
        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());
        
        tfEditor2D = new TransferFunction2DEditor(volume, gradients);
        tfEditor2D.addTFChangeListener(this);

        System.out.println("Finished initialization of RaycastRenderer");
    }

    public RaycastRendererPanel getPanel() {
        return panel;
    }

    public TransferFunction2DEditor getTF2DPanel() {
        return tfEditor2D;
    }
    
    public TransferFunctionEditor getTFPanel() {
        return tfEditor;
    }
     

    short getVoxel(double[] coord) {

        if (coord[0] < 0 || coord[0] > volume.getDimX() || coord[1] < 0 || coord[1] > volume.getDimY()
                || coord[2] < 0 || coord[2] > volume.getDimZ()) {
            return 0;
        }

        int x = (int) Math.floor(coord[0]);
        int y = (int) Math.floor(coord[1]);
        int z = (int) Math.floor(coord[2]);

        return volume.getVoxel(x, y, z);
    }
    
    //Tri-linear interpolation
    //Algorithm on slide 2-7
    //There is another algorithm mentioned at https://en.wikipedia.org/wiki/Trilinear_interpolation#Method
    short getVoxelTriLinear(double[] coord) {
        //Check validity
        //Since ceiling can be out of bound, so here we set to max - 1
        if (coord[0] < 0 || coord[0] > volume.getDimX() - 1 || coord[1] < 0 || coord[1] > volume.getDimY() - 1
                || coord[2] < 0 || coord[2] > volume.getDimZ() - 1) {
            return 0;
        }
        
        //Calculate x0, y0, z0 and x1, y1, z1
        int x0 = (int) Math.floor(coord[0]);
        int y0 = (int) Math.floor(coord[1]);
        int z0 = (int) Math.floor(coord[2]);
        
        int x1 = (int) Math.ceil(coord[0]);
        int y1 = (int) Math.ceil(coord[1]);
        int z1 = (int) Math.ceil(coord[2]);
        
        //Calculate the coefficiency alpha, beta and gamma
        double alpha = (coord[0] - x0) / (x1 - x0);
        double beta = (coord[1] - y0) / (y1 - y0);
        double gamma = (coord[2] - z0) / (z1 - z0);
        
        //We need to get the 8 cxxx points according to wiki
        //here 0 refers to x0, y0, z0 and 1 refers to x1, y1, z1
        //And c000 refers to SX0, c100 refers to SX1...c111 refers to SX7 in slide
        double c000 = volume.getVoxel(x0, y0, z0);
        double c001 = volume.getVoxel(x0, y0, z1);
        double c010 = volume.getVoxel(x0, y1, z0);
        double c011 = volume.getVoxel(x0, y1, z1);
        double c100 = volume.getVoxel(x1, y0, z0);
        double c101 = volume.getVoxel(x1, y0, z1);
        double c110 = volume.getVoxel(x1, y1, z0);
        double c111 = volume.getVoxel(x1, y1, z1);
        
        //Calculate from formula in slide
        
        return (short) Math.round(
                (1 - alpha) * (1 - beta) * (1 - gamma) * c000 +
                alpha * (1 - beta) * (1 - gamma) * c100 +
                (1 - alpha) * beta * (1 - gamma) * c010 +
                alpha * beta * (1 - gamma) * c110 +
                (1 - alpha) * (1 - beta) * gamma * c001 +
                alpha * (1 - beta) * gamma * c101 + 
                (1 - alpha) * beta * gamma * c011 + 
                alpha * beta * gamma * c111
        );  
    }
    
    void slicer(double[] viewMatrix) {

        // clear image
//        for (int j = 0; j < image.getHeight(); j++) {
//            for (int i = 0; i < image.getWidth(); i++) {
//                image.setRGB(i, j, 0);
//            }
//        }
        clearImage();

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

        
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2];

//                int val = getVoxel(pixelCoord); 
                //No longer use getVoxel but getVoxelTriLinear
                int val = getVoxelTriLinear(pixelCoord);
                
                // Map the intensity to a grey value by linear scaling
//                voxelColor.r = val/max;
//                voxelColor.g = voxelColor.r;
//                voxelColor.b = voxelColor.r;
//                voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                 voxelColor = tFunc.getColor(val);
                
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
//                image.setRGB(i, j, -0x10000000);
//                System.out.println("pixelcolor: " + pixelColor);
            }
        }

    }
    
    //MIP function
    //MIP is an extension of slicer
    //Slicer calculates the sum of St and MIP find tha max of it. 
    //So this method uses the structure of slicer
    //Another supplement method is added for calculating max st
    void MIP(double[] viewMatrix) {
        
        //clear image
        clearImage();
        
        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);
        
        // image is square
        int imageCenter = image.getWidth() / 2;
        int imageWidth = image.getWidth();
        int imageHeight = image.getHeight();
        int diagonal = (int) Math.ceil(
                    Math.sqrt((volume.getDimX() * volume.getDimX()) + (volume.getDimY() * volume.getDimY()) +
                                (volume.getDimZ() * volume.getDimZ())));

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);
        
        //Specific pixelCoord for slicer
        double[] pixelCoordSlicer = new double[3];
        
        // sample on a plane through the origin of the volume data
//        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();
        
        //We should also take mouse interaction into account
        //When pulling the image, it can be rendered in low resolution
        //When static, it can render in high resolution.
        //This can be changed by changing the MIPSlices and MIPSliceStep
        if (! interactiveMode) {
            MIPSliceStep = 1;
        } else {
            MIPSliceStep = 50;
        }

        
        for (int j = 1; j < imageHeight; j ++) {
            for (int i = 1; i < imageWidth; i ++) {
                //Max needed
                int maxVal = 0;
                //Basic slicer coord, center
                pixelCoordSlicer[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0];
                pixelCoordSlicer[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1];
                pixelCoordSlicer[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2];
                
                //calculate the other slices 
                //Here the max k should be limited to a number.
                //In principle, we should use either width or height for more slices.
                //However, it is too slow, here we hard code a number. 
                for (int k = 0; k < diagonal; k += MIPSliceStep) {
                    pixelCoord[0] = pixelCoordSlicer[0] + (k - (diagonal / 2)) * viewVec[0];
                    pixelCoord[1] = pixelCoordSlicer[1] + (k - (diagonal / 2)) * viewVec[1];
                    pixelCoord[2] = pixelCoordSlicer[2] + (k - (diagonal / 2)) * viewVec[2];
                    //old += does not work here, because the original pixelCoord will be added many times
                            
//                    int val = getVoxel(pixelCoord); 
                    //No longer use getVoxel but getVoxelTriLinear
                    int val = getVoxelTriLinear(pixelCoord);
                    //Check maxVal
                    if (val > maxVal) {
                        maxVal = val;
                    }
                }
//               
                
                // Map the intensity to a grey value by linear scaling
//                voxelColor.r = val/max;
//                voxelColor.r = maxVal / max;
//                voxelColor.g = voxelColor.r;
//                voxelColor.b = voxelColor.r;
//                voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
//                voxelColor.a = maxVal > 0 ? 1.0 : 0.0;
                // Alternatively, apply the transfer function to obtain a color
                // voxelColor = tFunc.getColor(val);
                //We use tFunc.getColor to get color on maxVal
                voxelColor = tFunc.getColor(maxVal);
                
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
//                image.setRGB(i, j, -0x10000000);
//                System.out.println("pixelcolor: " + pixelColor);
            }
        }
    }
    
    //Compositing
    void compositing(double[] viewMatrix) {
        clearImage();
        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;
        int imageWidth = image.getWidth();
        int imageHeight = image.getHeight();

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        int diagonal = (int) Math.ceil(
                    Math.sqrt((volume.getDimX() * volume.getDimX()) + (volume.getDimY() * volume.getDimY()) +
                                (volume.getDimZ() * volume.getDimZ())));
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);
        
        //Specific pixelCoord for slicer
        double[] pixelCoordSlicer = new double[3];

        // sample on a plane through the origin of the volume data
//        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();
        int step;
        
        if (! interactiveMode) {
            step = 1;
        } else {
            step = 50;
        }
        
        for (int j = 0; j < imageHeight; j++) {
            for (int i = 0; i < imageWidth; i++) {
                double r = 0;
                double g = 0;
                double b = 0;  
                pixelCoordSlicer[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                    + volumeCenter[0];
                pixelCoordSlicer[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                    + volumeCenter[1];
                pixelCoordSlicer[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                    + volumeCenter[2];
                
                for (int k = 0; k < diagonal; k += step) {
                    pixelCoord[0] = pixelCoordSlicer[0] + (k - (diagonal / 2)) * viewVec[0];
                    pixelCoord[1] = pixelCoordSlicer[1] + (k - (diagonal / 2)) * viewVec[1];
                    pixelCoord[2] = pixelCoordSlicer[2] + (k - (diagonal / 2)) * viewVec[2];

                    int val = getVoxelTriLinear(pixelCoord);
                    // Get color via compositing function
//                    voxelColor = new TFColor(0, 0, 0, 1);
                    TFColor bkgColor = tFunc.getColor(val);
                    r = (bkgColor.r * bkgColor.a) + ((1 - bkgColor.a) * r);
                    g = (bkgColor.g * bkgColor.a) + ((1 - bkgColor.a) * g);
                    b = (bkgColor.b * bkgColor.a) + ((1 - bkgColor.a) * b);

                    voxelColor = new TFColor(r, g, b, 1);
                }
                
                
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
//                image.setRGB(i, j, -0x10000000);
//                System.out.println("pixelcolor: " + pixelColor);
            }
        }
        
    }
    
    //2D Transfer function
    void transfer2D(double[] viewMatrix) {
        //clear image
        clearImage();
        
        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);
        
        // image is square
        int imageCenter = image.getWidth() / 2;
        int imageWidth = image.getWidth();
        int imageHeight = image.getHeight();
        
        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);
        
        TFColor voxelColor = new TFColor();
        
        //Specific pixelCoord for slicer
        double[] pixelCoordSlicer = new double[3];
        
        int diagonal = (int) Math.ceil(
                    Math.sqrt((volume.getDimX() * volume.getDimX()) + (volume.getDimY() * volume.getDimY()) +
                                (volume.getDimZ() * volume.getDimZ())));
        int step;
        
        //From the widget we have intensity, radius and opacity, as well as color
        //So we get the data first from the widget
        //Looks like TransferFunction2DEditor.trangleWidget.baseIntensity refers to indensity
        //radius refers to radius and color refers to color set
//        System.out.println("intensity: " + tfEditor2D.triangleWidget.baseIntensity);
        int setIntensity = tfEditor2D.triangleWidget.getIntensity();
        double radius = tfEditor2D.triangleWidget.getRadius();
        TFColor color = tfEditor2D.triangleWidget.getColor();
//        System.out.println("color: " + color.toString());

//        double opacity = 1;
        double voxelOpacity;
        
        if (! interactiveMode) {
            step = 1;
        } else {
            step = 50;
        }
        
        //We start with normal slicer
        for (int j = 0; j < imageHeight; j ++) {
            for (int i = 0 ; i < imageWidth; i ++) {
                double opacity = 1;
                voxelColor = color;
                //Basic slicer coord, center
                pixelCoordSlicer[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0];
                pixelCoordSlicer[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1];
                pixelCoordSlicer[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2];
                
                for (int k = 0; k < diagonal; k += step) {
                    pixelCoord[0] = pixelCoordSlicer[0] + (k - (diagonal / 2)) * viewVec[0];
                    pixelCoord[1] = pixelCoordSlicer[1] + (k - (diagonal / 2)) * viewVec[1];
                    pixelCoord[2] = pixelCoordSlicer[2] + (k - (diagonal / 2)) * viewVec[2];
                    
                    int val = getVoxelTriLinear(pixelCoord);
//                    int val = getVoxel2(pixelCoord);
                    
                    //Calculate 2D transfer
                    if (pixelCoord[0] < volume.getDimX() && pixelCoord[0] >= 0 &&
                            pixelCoord[1] < volume.getDimY() && pixelCoord[1] >= 0 &&
                            pixelCoord[2] < volume.getDimZ() && pixelCoord[2] >= 0) {
                        voxelOpacity = calculateOpacity(pixelCoord, val, setIntensity, radius, opacity);
                        if (shadeOn && voxelOpacity > 0) {
                            voxelColor = phoneShading(viewVec, pixelCoord, new TFColor(color.r, color.g, color.b, color.a));
                        }
                        //Apply the product function, first step
                        opacity *= 1 - voxelOpacity;
                    }
                }
                
                //Apply the rest part
                voxelColor.a = 1 - opacity;
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
//                image.setRGB(i, j, -0x10000000);
//                System.out.println("pixelcolor: " + pixelColor);
            }
        }
    }
    
    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glDisable(GL.GL_LINE_SMOOTH);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();

    }

    @Override
    public void visualize(GL2 gl) {


        if (volume == null) {
            return;
        }

        drawBoundingBox(gl);

        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);

        long startTime = System.currentTimeMillis();
        
        //Select main model
        if (this.status.equalsIgnoreCase("slicer")) {
            slicer(viewMatrix);
        } else if (this.status.equalsIgnoreCase("MIP")) {
            MIP(viewMatrix);
        } else if (this.status.equalsIgnoreCase("compositing")) {
            compositing(viewMatrix);
        } else if (this.status.equalsIgnoreCase("2DTrans")) {
            transfer2D(viewMatrix);
        }
                
//        slicer(viewMatrix);    
        
        long endTime = System.currentTimeMillis();
        double runningTime = (endTime - startTime);
        panel.setSpeedLabel(Double.toString(runningTime));

        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);

        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        // draw rendered image as a billboard texture
        texture.enable(gl);
        texture.bind(gl);
        double halfWidth = image.getWidth() / 2.0;
        gl.glPushMatrix();
        gl.glLoadIdentity();
        gl.glBegin(GL2.GL_QUADS);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glTexCoord2d(0.0, 0.0);
        gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
        gl.glTexCoord2d(0.0, 1.0);
        gl.glVertex3d(-halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 1.0);
        gl.glVertex3d(halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 0.0);
        gl.glVertex3d(halfWidth, -halfWidth, 0.0);
        gl.glEnd();
        texture.disable(gl);
        texture.destroy(gl);
        gl.glPopMatrix();

        gl.glPopAttrib();


        if (gl.glGetError() > 0) {
            System.out.println("some OpenGL error: " + gl.glGetError());
        }

    }
    private BufferedImage image;
    private double[] viewMatrix = new double[4 * 4];

    @Override
    public void changed() {
        for (int i=0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }
    
    //Supporting methods
    public void setStatus(String status) {
        this.status = status;
    }
    
    public void toggleShading() {
        this.shadeOn = !shadeOn;
//        System.out.println(shadeOn);
    }
    
    private void clearImage() {
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }
    }
    
    private double calculateOpacity(double[] pixelCoord, int voxelIntensity, int setIntensity, double radius, double alpha) {
        
        //Calculate gradientMagn
        VoxelGradient voxelGradient = gradients.getGradient((int) Math.floor(pixelCoord[0]), 
                (int) Math.floor(pixelCoord[1]), 
                (int) Math.floor(pixelCoord[2]));
        
        if (voxelGradient.mag < tfEditor2D.triangleWidget.getMinGradient() || 
                voxelGradient.mag > tfEditor2D.triangleWidget.getMaxGradient()) {
            return 0.0;
        }
        
        
        //Apply Levoy's equation 3
        //gradientMagn is mag of voxelGradient
        
        if (Math.abs(voxelGradient.mag) > gradients.getMaxGradientMagnitude()) {
            return 0.0;
        }
        
        if (Math.abs(voxelGradient.mag) == 0 && voxelIntensity == setIntensity) {
            return alpha;
        } else if (Math.abs(voxelGradient.mag) > 0 && 
                voxelIntensity - radius * Math.abs(voxelGradient.mag) <= setIntensity && 
                setIntensity <= voxelIntensity + radius * Math.abs(voxelGradient.mag)) {
            return alpha * (1 - (1 / radius * Math.abs((setIntensity - voxelIntensity) / 
                                                            Math.abs(voxelGradient.mag))));
        } else {
            return 0;
        }
    }
    
    // simplified phone shading with L = V
    //Ia = 0.1, kdiff = 0.7, kspec = 0.2, a = 10
    //I = Ia + Id * kdiff * (L * N) + kspec * (N * H)^a
    private TFColor phoneShading (double[] view, double[] pixelCoord, TFColor color) {
        //Calculate gradientMagn
        VoxelGradient voxelGradient = gradients.getGradient((int) Math.floor(pixelCoord[0]), 
                (int) Math.floor(pixelCoord[1]), 
                (int) Math.floor(pixelCoord[2]));
        //Normalize the gradient for N
        double[] N = new double[3];
        double gradientLength = Math.sqrt(voxelGradient.x * voxelGradient.x +
                voxelGradient.y * voxelGradient.y + voxelGradient.z * voxelGradient.z);
        VectorMath.setVector(N, voxelGradient.x / gradientLength,
                                    voxelGradient.y / gradientLength,
                                    voxelGradient.z / gradientLength);
        //Make the viewVec point away from the surface
        double[] viewVector = {-view[0], -view[1], -view[2]};
//        double viewVectorLength = VectorMath.length(viewVector);
        
        //L = V
        //L + V is two V
        double[] LPlusV = new double[3];
        VectorMath.setVector(LPlusV, viewVector[0] * 2, viewVector[1] * 2, viewVector[2] * 2);
        
        double LPlusVLength = VectorMath.length(LPlusV);

        //halfway vector
        double[] H = new double[3];
        H[0] = LPlusV[0] / LPlusVLength;
        H[1] = LPlusV[1] / LPlusVLength;
        H[2] = LPlusV[2] / LPlusVLength;
        
        double lDotN = VectorMath.dotproduct(viewVector, N);
        double hDotN = VectorMath.dotproduct(H, N);
        
        //Calculate new color
        TFColor newColor = new TFColor();

        if (lDotN <= 0) {
            return color; //In case the wrong direction
        }
        else {
                newColor.a = color.a;
        	double temp = Ia + kspec * Math.pow(hDotN, a);
        	newColor.r = color.r * kdiff * lDotN + temp;
        	newColor.g = color.g * kdiff * lDotN + temp;
        	newColor.b = color.b * kdiff * lDotN + temp;

        }
    	return newColor;
    }
    
}
