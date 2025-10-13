/*
 * This file is part of MiToBo, the Microscope Image Analysis Toolbox.
 *
 * Copyright (C) 2010 - @YEAR@
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Fore more information on MiToBo, visit
 *
 *    http://www.informatik.uni-halle.de/mitobo/
 *
 */

package de.unihalle.informatik.MiToBo.apps.plantCells.plastids;

import de.unihalle.informatik.Alida.exceptions.ALDOperatorException;
import de.unihalle.informatik.Alida.exceptions.ALDProcessingDAGException;
import de.unihalle.informatik.Alida.annotations.ALDAOperator;
import de.unihalle.informatik.Alida.annotations .ALDAOperator.Level;
import de.unihalle.informatik.Alida.annotations.Parameter;
import de.unihalle.informatik.Alida.datatypes.ALDFileString;
import de.unihalle.informatik.MiToBo.core.datatypes.MTBRegion2DSet;
import de.unihalle.informatik.MiToBo.core.datatypes.images.*;
import de.unihalle.informatik.MiToBo.core.datatypes.images.MTBImage.MTBImageType;
import de.unihalle.informatik.MiToBo.core.operator.*;
import de.unihalle.informatik.MiToBo.gui.MTBTableModel;
import de.unihalle.informatik.MiToBo.imageJ.plugins.cellCounter.datatypes.CellCntrMarker;
import de.unihalle.informatik.MiToBo.imageJ.plugins.cellCounter.datatypes.CellCntrMarkerShapeLine;
import de.unihalle.informatik.MiToBo.imageJ.plugins.cellCounter.datatypes.CellCntrMarkerVector;
import de.unihalle.informatik.MiToBo.imageJ.plugins.cellCounter.xml.ReadXML;
import de.unihalle.informatik.MiToBo.segmentation.regions.labeling.LabelComponentsSequential;
import de.unihalle.informatik.MiToBo.core.datatypes.MTBBorder2D;
import de.unihalle.informatik.MiToBo.visualization.drawing.DrawStringToImage;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.IOException;
import java.util.HashMap;
import java.util.TreeSet;
import java.util.Vector;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import loci.common.StatusEvent;
import loci.common.StatusListener;
import loci.common.StatusReporter;

/**
 * Operator to measure distances between plastids and nuclei per cell and to
 * evaluate stromules orientation with regard to cell nucleus.
 * 
 * @author Birgit Moeller
 */
@ALDAOperator(genericExecutionMode=ALDAOperator.ExecutionMode.ALL,
		level=Level.APPLICATION)
public class PlastidNucleiDistanceAnalyzer2D 
		extends MTBOperator implements StatusReporter {

	/**
	 * Identifier for outputs in verbose mode.
	 */
	private final static String opIdentifier = 
			"[PlastidNucleiDistanceAnalyzer2D] ";

	/**
	 * File with segmentation data.
	 */
	@Parameter(label = "Segmentation Data", required = true, dataIOOrder = 0,
		direction = Parameter.Direction.IN, description = "Segmentation data.")
	private ALDFileString inFileSegData = null;
	
	/**
	 * Input image stack.
	 */
	@Parameter(label = "Input Image Stack", required = true, dataIOOrder = 1,
		direction = Parameter.Direction.IN, description = "Input image.")
	private transient MTBImage inStack = null;

	/**
	 * Channel with boundary information.
	 */
	@Parameter( label= "Boundary Channel", 
		required = true, dataIOOrder = 2, direction = Parameter.Direction.IN, 
		description = "Boundary channel, starting with 1.")
	private int cellBoundaryChannel = 0;

	/**
	 * Type of cell nuclei.
	 */
	@Parameter( label= "Type of Nuclei Markers", 
		required = true, dataIOOrder = 3, direction = Parameter.Direction.IN, 
		description = "Nuclei marker type.")
	private int markerTypeNuclei = -1;

	/**
	 * Type of plastids without stromules.
	 */
	@Parameter( label= "Type of Plastid Markers", 
		required = true, dataIOOrder = 4, direction = Parameter.Direction.IN, 
		description = "Plastid marker type.")
	private int markerTypePlastid = -1;

	/**
	 * Type of plastids with stromules. 
	 */
	@Parameter( label= "Type of Stromuli Markers", 
		required = true, dataIOOrder = 5,	direction = Parameter.Direction.IN, 
		description = "Stromuli marker type.")
	private int markerTypePlastidWithStromuli = -1;

	/**
	 * Type of stromuli orientation markers.
	 */
	@Parameter( label= "Type of Stromuli Orientation Markers", 
		required = true, dataIOOrder = 6, direction = Parameter.Direction.IN, 
		description = "Stromuli orientation marker type.")
	private int markerTypeStromuliOrientation = -1;

	/**
	 * Result table with plastid-to-nuclei distances.
	 */
	@Parameter(label = "Result Table Distances", dataIOOrder = 0,
		direction = Parameter.Direction.OUT, 
		description = "Result table with plastid distance statistics.")
	private MTBTableModel resultDistanceTable = null;

	/**
	 * Result table with stromule orientation data.
	 */
	@Parameter(label = "Result Table Stromule Orientations", dataIOOrder = 1,
		direction = Parameter.Direction.OUT, 
		description = "Result table with stromule orientations.")
	private MTBTableModel resultStromulesOrientationTable = null;

	/**
	 * Debug image to visualize measurements.
	 */
	@Parameter(label = "Result Image", dataIOOrder = 2, required = true,
		direction = Parameter.Direction.OUT, description = "Output image.")
	private MTBImageRGB resultImage = null;
	
	/** 
	 * Vector of installed {@link StatusListener} objects.
	 */
	protected Vector<StatusListener> m_statusListeners;

	/**
	 * Default constructor.
	 * @throws ALDOperatorException	Thrown in case of operate failure.
	 */
	public PlastidNucleiDistanceAnalyzer2D() throws ALDOperatorException {
		this.m_statusListeners = new Vector<StatusListener>(1);
	}		
	 
	@Override
	protected void operate() 
			throws ALDOperatorException, ALDProcessingDAGException {

		this.resultImage = null;
		
		// post ImageJ status
		String msg = opIdentifier + "reading segmentation data file...";	
		this.notifyListeners(new StatusEvent(msg));

		ReadXML rxml;
		try {
			rxml = new ReadXML(this.inFileSegData.getFileName());
		} catch (SAXException e) {
			return;
		} catch (IOException e) {
			return;
		} catch (ParserConfigurationException e){
			return;
		}

		int width = this.inStack.getSizeX();
		int height = this.inStack.getSizeY();
		
		MTBImageByte boundaryImg = (MTBImageByte)this.inStack.getImagePart(0, 0, 
			this.cellBoundaryChannel-1, 0, 0,	width, height, 1, 1, 1).convertType(
				MTBImageType.MTB_BYTE, true);
		// invert image: boundaries in black
		for (int y=0;y<height;++y)
			for (int x=0;x<width;++x)
				if (boundaryImg.getValueInt(x, y) > 0)
					boundaryImg.putValueInt(x, y, 0);
				else
					boundaryImg.putValueInt(x, y, 255);
		
		// label cells
		LabelComponentsSequential labler = 
			new LabelComponentsSequential(boundaryImg, true);
		labler.runOp();
		MTBImage labelImg = labler.getLabelImage();
		MTBRegion2DSet rset = labler.getResultingRegions();
		
		// eliminate cells along image borders
		TreeSet<Integer> touchingLabels = new TreeSet<Integer>();
		for (int y=0;y<height;++y) {
			touchingLabels.add(new Integer(labelImg.getValueInt(0, y)));
			touchingLabels.add(new Integer(labelImg.getValueInt(width-1, y)));
		}
		for (int x=0;x<width;++x) {
			touchingLabels.add(new Integer(labelImg.getValueInt(x, 0)));
			touchingLabels.add(new Integer(labelImg.getValueInt(x, height-1)));
		}

		for (int y=0;y<height;++y) {
			for (int x=0;x<width;++x) {
				if (touchingLabels.contains(new Integer(labelImg.getValueInt(x, y))))
					labelImg.putValueInt(x, y, 0);
			}
		}
		
		this.resultImage = (MTBImageRGB)MTBImage.createMTBImage(
				width, height, 1, 1, 1, MTBImageType.MTB_RGB);
		int gval;
		for (int y=0;y<height;++y) {
			for (int x=0;x<width;++x) {
				if ((gval = labelImg.getValueInt(x, y)) == 0) {
					this.resultImage.putValueR(x, y, 200);
					this.resultImage.putValueG(x, y, 200);
					this.resultImage.putValueB(x, y, 200);					
				}
				else {
					this.resultImage.putValueR(x, y, gval);
					this.resultImage.putValueG(x, y, gval);
					this.resultImage.putValueB(x, y, gval);
				}
			}
		}
		
		// collect labels and sizes of regions
		int maxLabel =0;
		Integer val;
		TreeSet<Integer> labels = new TreeSet<>();
		for (int y=0;y<height;++y) {
			for (int x=0;x<width;++x) {
				val = new Integer(labelImg.getValueInt(x, y));
				if (val.intValue() != 0) {
					labels.add(val);
				}
				if (val.intValue() > maxLabel) {
					maxLabel = val.intValue();
				}
			}
		}
		int[] sizes = new int[maxLabel+1];
		for (int y=0;y<height;++y) {
			for (int x=0;x<width;++x) {
				if (labelImg.getValueInt(x, y) != 0)
					++sizes[labelImg.getValueInt(x, y)];
			}
		}
		
		CellCntrMarkerVector nucleiMarkers = null;
		CellCntrMarkerVector plastidMarkers = null;
		CellCntrMarkerVector stromuliMarkers = null;
		CellCntrMarkerVector stromuliOrientationMarkers = null;
		
		@SuppressWarnings("unchecked")
		Vector<CellCntrMarkerVector> mVectors = rxml.readMarkerData();
		
		for (CellCntrMarkerVector cv: mVectors) {
			if (   this.markerTypeNuclei != -1
					&& cv.getType() == this.markerTypeNuclei)
				nucleiMarkers = cv;
			else if (   this.markerTypePlastid != -1
						   && cv.getType() == this.markerTypePlastid)
				plastidMarkers = cv;
			else if (   this.markerTypePlastidWithStromuli != -1
					     && cv.getType() == this.markerTypePlastidWithStromuli)
				stromuliMarkers = cv;			
			else if (   this.markerTypeStromuliOrientation != -1
			     	   && cv.getType() == this.markerTypeStromuliOrientation)
				stromuliOrientationMarkers = cv;			
		}
		
		if (    nucleiMarkers == null 
				|| (plastidMarkers == null && stromuliMarkers == null))
			return;
		
		DrawStringToImage sDraw = new DrawStringToImage();

		// assign to each labeled cell the corresponding nuclei
		HashMap<Integer, CellCntrMarker> cellNuclei = new HashMap<>();
		Vector<Integer> multipleNucleiCells = new Vector<>();
		HashMap<Integer, Double> cellMaxDist = new HashMap<>();
		for (CellCntrMarker cm: nucleiMarkers) {
			int x = cm.getX();
			int y = cm.getY();
			int label = labelImg.getValueInt(x, y);
			// skip nuclei in the background
			if (label == 0)
				continue;
			this.resultImage.drawCircle2D(x, y, 0, 3, 0x00ff0000);
			if (!cellNuclei.containsKey(new Integer(label))) {
				cellNuclei.put(new Integer(label), cm);
				MTBBorder2D nmb = cm.getShape().getOutline();
				
				// calculate maximal possible distance in cell as maximum of 
				// minimal distances of all cell boundary points
				MTBBorder2D rb = rset.get(label-1).getBorder();
				double maxDist = 0, minDist;
				Point2D.Double maxc = null, maxn = null;
				for (Point2D.Double cp : rb.getPoints()) {
					Point2D.Double minc = null, minn = null;
					minDist = Double.MAX_VALUE;
					for (Point2D.Double np : nmb.getPoints()) {
						double dist = cp.distance(np);
						if (dist < minDist) {
							minDist = dist;
							minc = cp;
							minn = np;
						}
					}
					if (minDist > maxDist) {
						maxDist = minDist;
						maxc = minc;
						maxn = minn;
					}
				}
				cellMaxDist.put(new Integer(label), new Double(maxDist));
				this.resultImage.drawLine2D((int)maxc.x, (int)maxc.y, 
					(int)maxn.x, (int)maxn.y,	0x00aaaaaa);
				
				for (int i=1; i<nmb.getPointNum(); ++i) {
					Point2D.Double pp = nmb.getPointAt(i-1);
					Point2D.Double p = nmb.getPointAt(i);
					this.resultImage.drawLine2D(
							(int)pp.x, (int)pp.y, (int)p.x, (int)p.y, 0xff0000);
				}
				Point2D.Double pp = nmb.getPointAt(nmb.getPointNum()-1);
				Point2D.Double p = nmb.getPointAt(0);
				this.resultImage.drawLine2D(
						(int)pp.x, (int)pp.y, (int)p.x, (int)p.y, 0xff0000);

				sDraw.setInImage(this.resultImage);
				sDraw.setColor(Color.red);
				sDraw.setParameter("position", new Point2D.Double(x,y));
				sDraw.setString(Integer.toString(label));
				sDraw.runOp();
				this.resultImage = (MTBImageRGB)sDraw.getResultImage();
			}
			else
				multipleNucleiCells.add(new Integer(label));
		}
		
		int numPlastidsTotal = 0;
		int numStromuliTotal = 0;
		
		// iterate over all plastids and collect distance to nucleus
		HashMap<Integer, Vector<Double>> plastidDistPerCell = new HashMap<>();
		int x,y,label;
		Integer labelInt;
		if (plastidMarkers != null) {
			for (CellCntrMarker cm: plastidMarkers) {
				x = cm.getX();
				y = cm.getY();
				label = labelImg.getValueInt(x, y);
				labelInt = new Integer(label);

				// skip plastids in the background and in cells with multiple markers
				if (label == 0 || multipleNucleiCells.contains(labelInt))
					continue;
				
				CellCntrMarker nm = cellNuclei.get(new Integer(label));
				if (nm == null)
					continue;
				
				this.resultImage.drawCircle2D(x, y, 0, 1, 0x00ffff00);
				MTBBorder2D cmb = cm.getShape().getOutline();
				MTBBorder2D nmb = nm.getShape().getOutline();
				double minDist = Double.MAX_VALUE;
				Point2D.Double minpc = null, minpn = null;
				for (Point2D.Double cp : cmb.getPoints()) {
					for (Point2D.Double np : nmb.getPoints()) {
						double dist = cp.distance(np);
						if (dist < minDist) {
							minDist = dist;
							minpc = cp;
							minpn = np;
						}
					}
				}
				for (int i=1; i<cmb.getPointNum(); ++i) {
					Point2D.Double pp = cmb.getPointAt(i-1);
					Point2D.Double p = cmb.getPointAt(i);
					this.resultImage.drawLine2D(
							(int)pp.x, (int)pp.y, (int)p.x, (int)p.y, 0xffff00);
				}
				Point2D.Double pp = cmb.getPointAt(cmb.getPointNum()-1);
				Point2D.Double p = cmb.getPointAt(0);
				this.resultImage.drawLine2D(
						(int)pp.x, (int)pp.y, (int)p.x, (int)p.y, 0xffff00);

				this.resultImage.drawLine2D((int)minpc.x, (int)minpc.y, 
						(int)minpn.x, (int)minpn.y,	0x00ffff00);
				if (plastidDistPerCell.get(labelInt) == null)
					plastidDistPerCell.put(labelInt, new Vector<Double>());
				plastidDistPerCell.get(labelInt).add(new Double(minDist));
				
				sDraw.setInImage(this.resultImage);
				sDraw.setColor(Color.yellow);
				sDraw.setParameter("position", new Point2D.Double(x,y));
				sDraw.setString(
						Integer.toString(plastidDistPerCell.get(labelInt).size()));
				sDraw.runOp();
				this.resultImage = (MTBImageRGB)sDraw.getResultImage();
				++numPlastidsTotal;
			}
		}

		// remember IDs of plastids with stromuli
		HashMap<CellCntrMarker, Integer> plastidStromuleMarkerIDs = new HashMap<>();

		// iterate over all plastids with stromuli and collect distance to nucleus
		HashMap<Integer, Vector<Double>> stromuliDistPerCell = new HashMap<>();
		if (stromuliMarkers != null) {
			for (CellCntrMarker cm: stromuliMarkers) {
				x = cm.getX();
				y = cm.getY();
				label = labelImg.getValueInt(x, y);
				labelInt = new Integer(label);
				
				// skip plastids in the background
				if (label == 0 || multipleNucleiCells.contains(labelInt))
					continue;
				
				CellCntrMarker nm = cellNuclei.get(new Integer(label));
				if (nm == null)
					continue;
				
				this.resultImage.drawCircle2D(x, y, 0, 1, 0x0000ff00);
				MTBBorder2D cmb = cm.getShape().getOutline();
				MTBBorder2D nmb = nm.getShape().getOutline();
				double minDist = Double.MAX_VALUE;
				Point2D.Double minpc = null, minpn = null;
				for (Point2D.Double cp : cmb.getPoints()) {
					for (Point2D.Double np : nmb.getPoints()) {
						double dist = cp.distance(np);
						if (dist < minDist) {
							minDist = dist;
							minpc = cp;
							minpn = np;
						}
					}
				}
				for (int i=1; i<cmb.getPointNum(); ++i) {
					Point2D.Double pp = cmb.getPointAt(i-1);
					Point2D.Double p = cmb.getPointAt(i);
					this.resultImage.drawLine2D(
						(int)pp.x, (int)pp.y, (int)p.x, (int)p.y, 0x00ff00);
				}
				Point2D.Double pp = cmb.getPointAt(cmb.getPointNum()-1);
				Point2D.Double p = cmb.getPointAt(0);
				this.resultImage.drawLine2D(
						(int)pp.x, (int)pp.y, (int)p.x, (int)p.y, 0x00ff00);

				
				this.resultImage.drawLine2D((int)minpc.x, (int)minpc.y, 
						(int)minpn.x, (int)minpn.y,	0x0000ff00);
				if (stromuliDistPerCell.get(labelInt) == null)
					stromuliDistPerCell.put(labelInt, new Vector<Double>());
				stromuliDistPerCell.get(labelInt).add(new Double(minDist));
				
				sDraw.setInImage(this.resultImage);
				sDraw.setColor(Color.green);
				sDraw.setParameter("position", new Point2D.Double(x,y));
				sDraw.setString(
						Integer.toString(stromuliDistPerCell.get(labelInt).size()));
				sDraw.runOp();
				this.resultImage = (MTBImageRGB)sDraw.getResultImage();
				++numStromuliTotal;
				
				// remember the ID
				plastidStromuleMarkerIDs.put(cm, 
						new Integer(stromuliDistPerCell.get(labelInt).size()));
			}
		}

		
		// process stromuli orientation data
		if (stromuliOrientationMarkers != null && stromuliMarkers != null) {
			
			// init result table
			Vector<String> header = new Vector<>();
			header.add("Image-ID");
			header.add("Stromule-ID");
			header.add("Cell-ID");
			header.add("Plastid-ID"); 	
			header.add("Length"); 	
			header.add("Angle-diff");
			this.resultStromulesOrientationTable = new MTBTableModel(
					stromuliOrientationMarkers.size(), header.size(), header);

			int rowID = 0;
			Integer cellID;
			double angleDiff, orientAngle, nucleiAngle;
			for (CellCntrMarker cm: stromuliOrientationMarkers) {
				
				// check to cell the orientation belongs 
				CellCntrMarkerShapeLine cl = (CellCntrMarkerShapeLine)cm.getShape();
				double sx = cl.getStartPoint().x;
				double sy = cl.getStartPoint().y;
				double ex = cl.getEndPoint().x;
				double ey = cl.getEndPoint().y;
				label = labelImg.getValueInt((int)sx, (int)sy);
				labelInt = new Integer(label);
				cellID = labelInt;

				this.resultImage.drawLine2D(
					(int)sx, (int)sy, (int)ex, (int)ey, 0x0af2ef);
				this.resultImage.drawPoint2D((int)sx, (int)sy, 0, 0xff0000, 2);
				this.resultImage.drawPoint2D((int)ex, (int)ey, 0, 0x0000ff, 2);				

				// skip orientation markers in the background
				if (label == 0 || multipleNucleiCells.contains(labelInt))
					continue;
				
				CellCntrMarker nm = cellNuclei.get(new Integer(label));
				if (nm == null)
					continue;

				// get coordinates of nuclei center
				int ncx = nm.getX();
				int ncy = nm.getY();

				this.resultImage.drawLine2D(
						(int)sx, (int)sy, ncx, ncy, 0x0ab3f2);

				// check to which plastid the orientation marker belongs
				double minDist = Double.MAX_VALUE;
				CellCntrMarker minCP = null;
				for (CellCntrMarker cp: stromuliMarkers) {
					x = cp.getX();
					y = cp.getY();
					label = labelImg.getValueInt(x, y);
					labelInt = new Integer(label);

					// skip plastids in the background and in cells with multiple markers
					if (label == 0 || multipleNucleiCells.contains(labelInt))
						continue;

					CellCntrMarker cnm = cellNuclei.get(new Integer(label));
					if (cnm == null)
						continue;

					MTBBorder2D cpb = cp.getShape().getOutline();
					for (Point2D.Double p : cpb.getPoints()) {
						double dist = p.distance(cl.getStartPoint());
						if (dist < minDist) {
							minDist = dist;
							minCP = cp;
						}
					}
				}

				// angle between (start, nucleiCenter) and (start, end),
				// measured in degrees
				orientAngle = Math.toDegrees(Math.atan2((ey-sy),(ex-sx)));
				if (orientAngle < 0)
					orientAngle += 360;
				nucleiAngle = Math.toDegrees(Math.atan2((ncy-sy), (ncx-sx)));
				if (nucleiAngle < 0)
					nucleiAngle += 360;
				
				angleDiff = Math.abs(orientAngle-nucleiAngle);
				if (angleDiff > 180)
					angleDiff -= 180;
							
				// image name
				this.resultStromulesOrientationTable.setValueAt(
					this.inStack.getTitle(), rowID, 0);
				// stromule ID
				this.resultStromulesOrientationTable.setValueAt(
					new Integer(rowID+1), rowID, 1);	
				// cell ID to which the stromule belongs
				this.resultStromulesOrientationTable.setValueAt(cellID, rowID, 2);	
				// plastid ID to which the stromule belongs
				this.resultStromulesOrientationTable.setValueAt(
					plastidStromuleMarkerIDs.get(minCP), rowID, 3);
				// estimated length
				this.resultStromulesOrientationTable.setValueAt(
					new Double(cl.getStartPoint().distance(cl.getEndPoint())), rowID, 4);
				// orientation difference
				this.resultStromulesOrientationTable.setValueAt(
					new Double(angleDiff), rowID, 5);
				++rowID;		

			} // end of for-loop over all stromule orientation markers
 		}

		Vector<String> header = new Vector<>();
		header.add("Image-ID");
		header.add("Cell-ID");
		header.add("Plastid-ID"); 	
		header.add("HasStromule");
		header.add("Distance");
		header.add("Cell-Area"); 	
		header.add("Summ-plastids");
		header.add("Max-pos-dist");
		this.resultDistanceTable = new MTBTableModel(
				numPlastidsTotal+numStromuliTotal, header.size(), header);
		
		int rowID = 0;
		// iterate over all existing labels
		for(Integer l : labels) {
			
			// cells without nucleus
			if (cellNuclei.get(l) == null) {
				this.resultDistanceTable.setValueAt(this.inStack.getTitle(), rowID, 0);
				this.resultDistanceTable.setValueAt(l, rowID, 1);				
				this.resultDistanceTable.setValueAt("no nucleus found", rowID, 2);
				this.resultDistanceTable.setValueAt("undefined", rowID, 3);
				this.resultDistanceTable.setValueAt("undefined", rowID, 4);
				this.resultDistanceTable.setValueAt(
						new Integer(sizes[l.intValue()]), rowID, 5);
				this.resultDistanceTable.setValueAt("undefined", rowID, 6);
				this.resultDistanceTable.setValueAt("undefined", rowID, 7);
				++rowID;		
				continue;
			}
			
			// cells with multiple nuclei
			if (multipleNucleiCells.contains(l)) {
				this.resultDistanceTable.setValueAt(this.inStack.getTitle(), rowID, 0);
				this.resultDistanceTable.setValueAt(l, rowID, 1);				
				this.resultDistanceTable.setValueAt(
						"more than one nucleus found", rowID, 2);
				this.resultDistanceTable.setValueAt("undefined", rowID, 3);
				this.resultDistanceTable.setValueAt("undefined", rowID, 4);
				this.resultDistanceTable.setValueAt(
						new Integer(sizes[l.intValue()]), rowID, 5);
				this.resultDistanceTable.setValueAt("undefined", rowID, 6);
				this.resultDistanceTable.setValueAt("undefined", rowID, 7);
				++rowID;
				continue;
			}
			
			// cells with a single nucleus
			int pID = 1;

			Vector<Double> pdists = plastidDistPerCell.get(l);
			int numPlastids = 0;
			if (pdists != null) 
				numPlastids = pdists.size();
			
			Vector<Double> sdists = stromuliDistPerCell.get(l);
			int numStromuli = 0;
			if (sdists != null) 
				numStromuli = sdists.size();
			
			// total number of detected plastids with and without stromuli in cell l
			int numPlastidsAndStromuli = numPlastids + numStromuli;
			
			if (pdists != null) {
				for (Double d: pdists) {
					this.resultDistanceTable.setValueAt(this.inStack.getTitle(), rowID, 0);
					this.resultDistanceTable.setValueAt(l, rowID, 1);				
					this.resultDistanceTable.setValueAt(new Integer(pID), rowID, 2);
					this.resultDistanceTable.setValueAt(new Integer(0), rowID, 3);
					this.resultDistanceTable.setValueAt(d, rowID, 4);
					this.resultDistanceTable.setValueAt(
							new Integer(sizes[l.intValue()]), rowID, 5);
					this.resultDistanceTable.setValueAt(
							new Integer(numPlastidsAndStromuli), rowID, 6);
					this.resultDistanceTable.setValueAt(cellMaxDist.get(l), rowID, 7);
					++pID;
					++rowID;
				}
			}
			if (sdists != null) {
				for (Double d: sdists) {
					this.resultDistanceTable.setValueAt(this.inStack.getTitle(), rowID, 0);
					this.resultDistanceTable.setValueAt(l, rowID, 1);				
					this.resultDistanceTable.setValueAt(new Integer(pID), rowID, 2);
					this.resultDistanceTable.setValueAt(new Integer(1), rowID, 3);
					this.resultDistanceTable.setValueAt(d, rowID, 4);
					this.resultDistanceTable.setValueAt(
							new Integer(sizes[l.intValue()]), rowID, 5);
					this.resultDistanceTable.setValueAt(
							new Integer(numPlastidsAndStromuli), rowID, 6);
					this.resultDistanceTable.setValueAt(cellMaxDist.get(l), rowID, 7);
					++pID;
					++rowID;
				}
			}
			
			if (pdists == null && sdists == null) {
				this.resultDistanceTable.setValueAt(this.inStack.getTitle(), rowID, 0);
				this.resultDistanceTable.setValueAt(l, rowID, 1);				
				this.resultDistanceTable.setValueAt(new Integer(0), rowID, 2);
				this.resultDistanceTable.setValueAt(new Integer(0), rowID, 3);
				this.resultDistanceTable.setValueAt("undefined", rowID, 4);
				this.resultDistanceTable.setValueAt(
						new Integer(sizes[l.intValue()]), rowID, 5);
				this.resultDistanceTable.setValueAt(new Integer(0), rowID, 6);
				this.resultDistanceTable.setValueAt(cellMaxDist.get(l), rowID, 7);
				++pID;
				++rowID;				
			}
		}
		this.resultImage.setTitle(
			"Plastid distance and stromule orientation statistics");
			
		msg = opIdentifier + "calcuations done!";	
		this.notifyListeners(new StatusEvent(msg));
	}
		
	// ----- StatusReporter interface
	
	@Override
	public void addStatusListener(StatusListener statuslistener) {	
		this.m_statusListeners.add(statuslistener);	
	}

	@Override
	public void notifyListeners(StatusEvent e) {
		for (int i = 0; i < this.m_statusListeners.size(); i++) {
			this.m_statusListeners.get(i).statusUpdated(e);
		}
	}

	@Override
	public void removeStatusListener(StatusListener statuslistener) {
		this.m_statusListeners.remove(statuslistener);
	}	
}

