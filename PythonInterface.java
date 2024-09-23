package Filtration;

import py4j.GatewayServer;

import topcat.matrix.distancematrix.DistanceMatrix;
import topcat.util.BinomialCoeffTable;
import topcat.util.IntTuple;
import topcat.util.Point;
import topcat.util.Grid;
import topcat.util.GridIterator;
import topcat.persistence.simplex.Simplex;
import topcat.persistence.simplex.SimplexStorageStructure;
import topcat.persistence.simplex.SimplicialComplex;

import it.unimi.dsi.fastutil.Arrays;
import it.unimi.dsi.fastutil.Function;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;

import javax.sound.midi.SysexMessage;
import javax.swing.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;

import java.util.Collections;
import java.util.Arrays;
import java.util.Locale;
import java.util.Scanner;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class Filtration {
    OwnSimplexStorageStructure storageStructure;
    String type;
    int numVertices;
    int max_dim;
    int length_hom;
    int length_zigzag;

    // List<DistanceMatrix>: Liste zur Zeit t_1, t_12, t_2, t_23, usw (t1, t2 usw
    // sogar unnötig, auslesen von blockmatrizen)
    // direkt in Datei schreiben, nicht als Liste speichern
    //
    public Filtration(List<List<Point>> timeseries, List<Double> filtrationVals, int max_dim) {
        List<DistanceMatrix> distanceMatrices = new ArrayList<>();
        if (timeseries.size() == 1) {
            DistanceMatrix distmat = DistanceMatrix.computeEuclideanDistanceMatrix(timeseries.get(0));
            distanceMatrices.add(distmat);
        }
        for (int i = 0; i < timeseries.size() - 1; i++) {
            List<Point> tmp_points = new ArrayList<>();
            tmp_points.addAll(timeseries.get(i));
            tmp_points.addAll(timeseries.get(i + 1));
            DistanceMatrix tmp = DistanceMatrix.computeEuclideanDistanceMatrix(tmp_points);
            distanceMatrices.add(tmp);
        }
        String type = "f";
        for (int i = 0; i < timeseries.size() - 1; i++) {
            type = type.concat("gf");
        }
        this.type = type;
        storageStructure = SimplicialComplexTimeseries.computeSimplexStream(timeseries, distanceMatrices,
                filtrationVals, max_dim);

        int num_vertices = 0;
        for (List<Point> window : timeseries) {
            num_vertices += window.size();
        }
        numVertices = num_vertices;
        length_zigzag = timeseries.size() * 2 - 1;
        length_hom = filtrationVals.size();
        this.max_dim = max_dim;

        System.out.println("Computing filtration completed!");

        // length_hom = storageStructure.gridSizeCopy.get(0)+1;
        // length_zigzag = storageStructure.gridSizeCopy.get(1);
    }

    public void writeToFiles(String filename) throws Exception {
        // String entire_filename = "/HOME/flammer/projects/fzz-main/build/data/";
        String entire_filename = "/LocalData/flammer/data_landscapes/";
        entire_filename = entire_filename.concat(filename).concat("_");
        // String filename = "test_";
        for (int i = 0; i < length_hom; i++) {
            for (int j = 0; j < length_zigzag; j++) {
                // for(int i=0; i<2; i++){
                // for (int j=0; j<2; j++){
                IntTuple x = new IntTuple(i, j);
                FiltrationAtX tmp = new FiltrationAtX(storageStructure, max_dim, length_hom, length_zigzag, x, type,
                        false);
                String tuple = Integer.toString(i).concat("_").concat(Integer.toString(j));
                tmp.writeSimplexwiseFiltration(entire_filename.concat(tuple).concat("_input").concat(".txt"), max_dim,
                        numVertices);
                tmp.simplexwiseToBarcode(entire_filename.concat(tuple).concat("_barcode_transl.txt"));
            }
        }
    }

    public void writeToFilesHomogeneousZigzag(String filename, IntTuple x) throws Exception {
        // String entire_filename = "/HOME/flammer/projects/fzz-main/build/data/";
        String entire_filename = "/LocalData/flammer/data_landscapes/";
        entire_filename = entire_filename.concat(filename).concat("_");
        FiltrationAtX tmp = new FiltrationAtX(storageStructure, max_dim, length_hom, length_zigzag, x, type, true);
        String tuple = Integer.toString(x.get(0)).concat("_").concat(Integer.toString(x.get(1)).concat("_homogeneous"));
        tmp.writeSimplexwiseFiltration(entire_filename.concat(tuple).concat("_input").concat(".txt"), max_dim,
                numVertices);
        tmp.simplexwiseToBarcode(entire_filename.concat(tuple).concat("_barcode_transl.txt"));
    }

    public class FiltrationAtX {
        HashMap<Integer, List<List<Simplex>>> filtrationAlongPath;
        List<Integer> relevantPoints;
        List<Character> insertion_deletion;

        FiltrationAtX(OwnSimplexStorageStructure storageStructure, int max_dim, int length_hom, int length_zigzag,
                IntTuple x, String type, boolean homogeneous_zigzag) throws Exception {
            List<Integer> int_list = new ArrayList<>();
            if (x.get(0) < 0 || x.get(1) < 0 || x.get(0) >= length_hom || x.get(1) >= length_zigzag) {
                throw new Exception("Intervals at index are not defined!");
            }
            int_list.add(x.get(0));
            int_list.add(x.get(1));
            int_list.add(length_hom - x.get(0) - 1);
            int_list.add(length_zigzag - x.get(1) - 1);
            IntTuple comp = new IntTuple(int_list);
            int radius_max = comp.min();

            List<IntTuple> path = new ArrayList<>();
            List<IntTuple> pointsOfInterest = new ArrayList<>();
            if (homogeneous_zigzag) {
                for (int z = 0; z < length_zigzag; z++) {
                    List<Integer> index = new ArrayList<>();
                    index.add(x.get(0));
                    index.add(z);
                    IntTuple ind = new IntTuple(index);
                    path.add(ind);
                }
                pointsOfInterest.addAll(path);
            } else {
                path.addAll(Path.getEntirePath(x, radius_max, true));
                System.out.println("entire path size:" + path.size());
                pointsOfInterest.addAll(Path.getPointsOfInterest(path, x));
            }

            // test:
            for (IntTuple point : pointsOfInterest) {
                System.out.println("Point of interest: " + point.toString());
            }

            List<Integer> relevant_points = new ArrayList<>();
            for (int i = 0; i < pointsOfInterest.size(); i++) {
                for (int j = 0; j < path.size(); j++) {
                    if (pointsOfInterest.get(i).equals(path.get(j))) {
                        relevant_points.add(j);
                        break;
                    }
                }
            }
            this.relevantPoints = relevant_points;
            String type_alon_path = Path.typeAlongPath(path, type);
            System.out.println("Type along path: " + type_alon_path);
            List<Character> insertion_deletion = new ArrayList<>();
            insertion_deletion.add('i');
            HashMap<Integer, List<List<Simplex>>> filtr = new HashMap<>();

            // initialisation:
            System.out.println("Initialisation start");
            for (int dim = 0; dim <= max_dim; dim++) {
                List<List<Simplex>> tmplist = new ArrayList<>();
                tmplist.add(storageStructure.getSimplicesLEQThan(dim, path.get(0)));
                filtr.put(dim, tmplist);
            }
            System.out.println("Initialisation end");

            for (int dim = 0; dim <= max_dim; dim++) {
                List<List<Simplex>> changesAlongPath = filtr.get(dim);
                for (int i = 0; i < path.size() - 1; i++) {
                    IntTuple from = path.get(i);
                    IntTuple to = path.get(i + 1);
                    List<Simplex> simplices1 = storageStructure.getSimplicesLEQThan(dim, from);
                    List<Simplex> simplices2 = storageStructure.getSimplicesLEQThan(dim, to);
                    List<Simplex> difference = OwnSimplexStorageStructure.differenceOfSimplices(simplices1, simplices2);
                    changesAlongPath.add(difference);
                }
                filtr.put(dim, changesAlongPath);
            }
            for (int i = 0; i < type_alon_path.length(); i++) {
                if (type_alon_path.charAt(i) == 'f') {
                    insertion_deletion.add('i');
                } else if (type_alon_path.charAt(i) == 'g') {
                    insertion_deletion.add('d');
                }
            }
            // last step:
            insertion_deletion.add('d');
            for (int dim = 0; dim < max_dim + 1; dim++) {
                List<List<Simplex>> tmplist = new ArrayList<>();
                tmplist.add(storageStructure.getSimplicesLEQThan(dim, path.get(path.size() - 1)));
                filtr.get(dim).addAll(tmplist);
            }

            System.out.println("Insertion deletion " + insertion_deletion);

            filtrationAlongPath = filtr;
            this.insertion_deletion = insertion_deletion;
            System.out.println("insertion deletion size: " + insertion_deletion.size());
        }

        public List<Integer> simplexwiseToBarcode(String filename) {
            List<Integer> out = new ArrayList<>();
            out.add(0);
            List<Integer> path_number_simplices = new ArrayList<>();
            // for(int i=0; i<insertion_deletion.size(); i++) {
            if (relevantPoints.size() > 1) {
                for (int i = 0; i < filtrationAlongPath.get(0).size(); i++) {
                    int sum = 0;
                    for (int dim = 0; dim < filtrationAlongPath.keySet().size(); dim++) {
                        sum += filtrationAlongPath.get(dim).get(i).size();
                    }
                    path_number_simplices.add(sum);
                }
                for (int j = 0; j < relevantPoints.size() - 1; j++) {
                    int start_incl = relevantPoints.get(j);
                    int end_excl = relevantPoints.get(j + 1);
                    int sum_new = 0;
                    for (int i = start_incl; i < end_excl; i++) {
                        sum_new += path_number_simplices.get(i);
                    }
                    sum_new += out.get(out.size() - 1);
                    out.add(sum_new);
                }
            }
            System.out.println("relevant points size: " + relevantPoints.size());
            for (Integer point : relevantPoints) {
                System.out.println(point);
            }

            // write List to file
            try {
                FileWriter myWriter = new FileWriter(filename);
                for (Integer i : out) {
                    myWriter.write(Integer.toString(i));
                    myWriter.write("\n");
                }
                myWriter.close();
                System.out.println("Successfully wrote to the file.");
            } catch (IOException e) {
                System.out.println("An error occurred.");
                e.printStackTrace();
            }
            return out;
        }

        public List<Simplex> writeSimplexwiseFiltration(String filename, int maxDim, int numVertices) {
            List<Simplex> out = new ArrayList<>();

            BinomialCoeffTable binomialCoeffTable = new BinomialCoeffTable(numVertices, maxDim);
            int sum = 0; // nur zum testen
            // String entire_filename = "/HOME/flammer/projects/fzz-main/build/data/";
            // entire_filename = entire_filename.concat(filename);
            try {
                FileWriter myWriter = new FileWriter(filename);

                for (int i = 0; i < insertion_deletion.size(); i++) {
                    if (insertion_deletion.get(i) == 'i') {
                        for (int dim = 0; dim <= maxDim; dim++) {
                            System.out.println("Filtration along path: " + filtrationAlongPath.get(dim).get(i).size()
                                    + ", dim: " + dim);
                            List<Simplex> list = filtrationAlongPath.get(dim).get(i);
                            for (Simplex simplex : list) {
                                int[] vertices = Simplex.get_simplex_vertices(simplex.getIndex(),
                                        simplex.getDimension(), numVertices, binomialCoeffTable);
                                myWriter.write(insertion_deletion.get(i).toString());
                                for (int k = 0; k < vertices.length; k++) {
                                    myWriter.write(" ");
                                    myWriter.write(Integer.toString(vertices[k]));

                                }
                                sum++;
                                myWriter.write("\n");
                            }
                        }
                    } else if (insertion_deletion.get(i) == 'd') {
                        for (int dim = maxDim; dim >= 0; dim--) {
                            List<Simplex> list = filtrationAlongPath.get(dim).get(i);
                            for (Simplex simplex : list) {
                                int[] vertices = Simplex.get_simplex_vertices(simplex.getIndex(),
                                        simplex.getDimension(), numVertices, binomialCoeffTable);
                                myWriter.write(insertion_deletion.get(i).toString());
                                for (int k = 0; k < vertices.length; k++) {
                                    myWriter.write(" ");
                                    myWriter.write(Integer.toString(vertices[k]));

                                }
                                sum++;
                                myWriter.write("\n");
                            }
                        }
                    }

                }
                myWriter.close();
                System.out.println("Successfully wrote to the file.");
            } catch (IOException e) {
                System.out.println("An error occurred.");
                e.printStackTrace();
            }

            System.out.println("sum simplices : " + sum);
            int sum2 = 0;
            for (int dim = 0; dim <= maxDim; dim++) {
                System.out.println("filtragtion along path: " + filtrationAlongPath.get(dim).size());
                int tmp_sum = 0;
                List<List<Simplex>> tmp = filtrationAlongPath.get(dim);
                for (List<Simplex> list : tmp) {
                    for (Simplex simplex : list) {
                        tmp_sum++;
                    }

                }
                sum2 += tmp_sum;
            }
            System.out.println("Sum 2: " + sum2);

            return out;
        }
    }

}

public class OwnSimplexStorageStructure extends SimplexStorageStructure {
    IntTuple gridSizeCopy;
    Int2ObjectOpenHashMap<Grid<List<Simplex>>> simplexContainerCopy;

    public OwnSimplexStorageStructure(List<List<Double>> filtrationValues, IntTuple gridSize, Integer max_dimesion,
            Integer n_vertices) {
        super(filtrationValues, gridSize, max_dimesion, n_vertices);
        this.gridSizeCopy = gridSize;
        simplexContainerCopy = new Int2ObjectOpenHashMap<>();
    }

    @Override
    public List<Simplex> getSimplicesLEQThan(int dim, IntTuple filtrationIndex) {
        List<Simplex> simplices = new ArrayList<>();
        List<IntTuple> indicesList = new ArrayList<>();
        if (filtrationIndex.get(1) % 2 == 0) {
            for (int i = 0; i < filtrationIndex.get(0) + 1; i++) {
                IntTuple index = new IntTuple(i, filtrationIndex.get(1));
                indicesList.add(index);
            }
        } else {
            for (int k = filtrationIndex.get(1) - 1; k <= filtrationIndex.get(1) + 1; k++) {
                for (int i = 0; i < filtrationIndex.get(0) + 1; i++) {
                    IntTuple index = new IntTuple(i, k);
                    indicesList.add(index);
                }
            }
        }

        for (IntTuple v : indicesList) {
            List<Simplex> local_simplices = getSimplicesAt(dim, v);
            // System.out.println("IntTuple " + v.toString() + " dim " + dim);
            if (local_simplices != null) {
                simplices.addAll(local_simplices);
            }
        }
        // System.out.println("Index list: " + indicesList);
        Collections.sort(simplices);
        return simplices;
    }

    @Override
    public void addElement(Simplex simplex, IntTuple filtrationIndex) {
        if (!simplexContainerCopy.containsKey(simplex.getDimension())) {
            simplexContainerCopy.put(simplex.getDimension(), Grid.create(gridSizeCopy));
        }
        List<Simplex> simplices = simplexContainerCopy.get(simplex.getDimension()).get(filtrationIndex);
        if (simplices == null) {
            simplices = new ArrayList<>();
            simplexContainerCopy.get(simplex.getDimension()).set(filtrationIndex, simplices);
        }
        simplices.add(simplex);
        // simplexContainerCopy.get(simplex.getDimension()).get(filtrationIndex).add(simplex);
    }

    // das sind nur die Simplices, die neu auftauchen?
    @Override
    public List<Simplex> getSimplicesAt(int dim, IntTuple filtrationIndex) {
        if (!simplexContainerCopy.containsKey(dim)) {
            return null;
        }
        return simplexContainerCopy.get(dim).get(filtrationIndex);
    }

    public List<Simplex> addSimplices(IntTuple from, IntTuple to, int dim) {
        List<Simplex> simplices = new ArrayList<>();
        List<Simplex> simplicesTo = getSimplicesAt(dim, to);
        List<Simplex> simplicesFrom = getSimplicesAt(dim, from);
        for (Simplex simplexTo : simplicesTo) {
            boolean equal = false;
            for (Simplex simplexFrom : simplicesFrom) {
                if (simplexFrom.equals(simplexTo)) {
                    equal = true;
                    break;
                }
            }
            if (!equal) {
                simplices.add(simplexTo);
            }
        }
        return simplices;
    }

    public List<Simplex> removeSimplices(IntTuple from, IntTuple to, int dim) {
        List<Simplex> simplices = new ArrayList<>();
        List<Simplex> simplicesTo = getSimplicesAt(dim, to);
        List<Simplex> simplicesFrom = getSimplicesAt(dim, from);
        for (Simplex simplexFrom : simplicesFrom) {
            boolean equal = false;
            for (Simplex simplexTo : simplicesTo) {
                if (simplexFrom.equals(simplexTo)) {
                    equal = true;
                    break;
                }
            }
            if (!equal) {
                simplices.add(simplexFrom);
            }
        }
        return simplices;
    }

    public static List<Simplex> differenceOfSimplices(List<Simplex> list1, List<Simplex> list2) {
        List<Simplex> out = new ArrayList<>();
        List<Simplex> smaller = new ArrayList<>();
        List<Simplex> bigger = new ArrayList<>();
        if (list1.size() < list2.size()) {
            smaller.addAll(list1);
            bigger.addAll(list2);
        } else {
            smaller.addAll(list2);
            bigger.addAll(list1);
        }
        for (Simplex simplex2 : bigger) {
            boolean equal = false;
            for (Simplex simplex1 : smaller) {
                if (simplex1.equals(simplex2)) {
                    equal = true;
                    break;
                }
            }
            if (!equal) {
                out.add(simplex2);
            }
        }
        return out;
    }
}

public class Path {
    public static List<IntTuple> getFirstPathIndices(IntTuple x, int radius, boolean totheright) {
        if (radius > 0 && totheright) {
            List<IntTuple> pathIndices = new ArrayList<>();
            IntTuple e1 = new IntTuple(Arrays.asList(0, 1));
            IntTuple ez = new IntTuple(Arrays.asList(1, 0));
            IntTuple tmp_index = new IntTuple(x);
            for (int i = 0; i < radius; i++) {
                tmp_index = tmp_index.minus(e1).minus(ez);
            }
            pathIndices.add(tmp_index);

            for (int i = 0; i < 2 * radius; i++) {
                tmp_index = tmp_index.plus(e1);
                pathIndices.add(tmp_index);
            }
            tmp_index = tmp_index.plus(ez);
            pathIndices.add(tmp_index);
            // tmp_index = tmp_index.minus(e1);
            // pathIndices.add(tmp_index);
            pathIndices.addAll(getFirstPathIndices(x, radius - 1, false));
            return pathIndices;
        } else if (radius > 0) {
            List<IntTuple> pathIndices = new ArrayList<>();
            IntTuple e1 = new IntTuple(Arrays.asList(0, 1));
            IntTuple ez = new IntTuple(Arrays.asList(1, 0));
            IntTuple tmp_index = new IntTuple(x);
            for (int i = 0; i < radius; i++) {
                tmp_index = tmp_index.plus(e1).minus(ez);
            }
            pathIndices.add(tmp_index);

            for (int i = 0; i < 2 * radius; i++) {
                tmp_index = tmp_index.minus(e1);
                pathIndices.add(tmp_index);
            }
            tmp_index = tmp_index.plus(ez);
            pathIndices.add(tmp_index);
            // tmp_index = tmp_index.minus(e1);
            // pathIndices.add(tmp_index);
            pathIndices.addAll(getFirstPathIndices(x, radius - 1, true));
            return pathIndices;
        } else {
            List<IntTuple> pathIndices = new ArrayList<>();
            IntTuple center = new IntTuple(Arrays.asList(x.get(0), x.get(1)));
            pathIndices.add(center);
            return pathIndices;
        }
    }

    public static List<IntTuple> getSecondPathIndices(IntTuple x, int radius, boolean totheleft) {
        if (radius > 0 && totheleft) {
            List<IntTuple> pathIndices = new ArrayList<>();
            IntTuple e1 = new IntTuple(Arrays.asList(0, 1));
            IntTuple ez = new IntTuple(Arrays.asList(1, 0));
            IntTuple tmp_index = new IntTuple(x);
            for (int i = 0; i < radius; i++) {
                tmp_index = tmp_index.plus(e1).plus(ez);
            }
            pathIndices.add(tmp_index);

            for (int i = 0; i < 2 * radius; i++) {
                tmp_index = tmp_index.minus(e1);
                pathIndices.add(tmp_index);
            }
            tmp_index = tmp_index.minus(ez);
            pathIndices.add(tmp_index);
            // tmp_index = tmp_index.minus(e1);
            // pathIndices.add(tmp_index);
            pathIndices.addAll(getSecondPathIndices(x, radius - 1, false));
            return pathIndices;
        } else if (radius > 0) {
            List<IntTuple> pathIndices = new ArrayList<>();
            IntTuple e1 = new IntTuple(Arrays.asList(0, 1));
            IntTuple ez = new IntTuple(Arrays.asList(1, 0));
            IntTuple tmp_index = new IntTuple(x);
            for (int i = 0; i < radius; i++) {
                tmp_index = tmp_index.minus(e1).plus(ez);
            }
            pathIndices.add(tmp_index);

            for (int i = 0; i < 2 * radius; i++) {
                tmp_index = tmp_index.plus(e1);
                pathIndices.add(tmp_index);
            }
            tmp_index = tmp_index.minus(ez);
            pathIndices.add(tmp_index);
            // tmp_index = tmp_index.minus(e1);
            // pathIndices.add(tmp_index);
            pathIndices.addAll(getSecondPathIndices(x, radius - 1, true));
            return pathIndices;
        } else {
            // pathIndices.add(x);
            return new ArrayList<>();
        }
    }

    public static List<IntTuple> reverseList(List<IntTuple> list) {
        List<IntTuple> newlist = new ArrayList<>();
        while (list.size() > 0) {
            newlist.add(list.get(list.size() - 1));
            list.remove(list.size() - 1);
        }
        return newlist;
    }

    public static List<IntTuple> getEntirePath(IntTuple x, int radius, boolean toTheRight) {
        List<IntTuple> out = new ArrayList<>();
        if (toTheRight) {
            out.addAll(getFirstPathIndices(x, radius, true));
            out.addAll(reverseList(getSecondPathIndices(x, radius, true)));
        } else {
            out.addAll(getFirstPathIndices(x, radius, false));
            out.addAll(reverseList(getSecondPathIndices(x, radius, false)));
        }
        return out;
    }

    public static List<IntTuple> getPointsOfInterest(List<IntTuple> indexlist, IntTuple x) {
        List<IntTuple> points = new ArrayList<>();
        points.add(indexlist.get(0));
        if (indexlist.size() == 1) {
            return points;
        }
        for (int i = 0; i < indexlist.size() - 1; i++) {
            if (indexlist.get(i).get(0) < indexlist.get(i + 1).get(0) && i + 2 < indexlist.size()
                    && indexlist.get(i).get(0) < x.get(0)) {
                points.add(indexlist.get(i + 2));
            } else if (indexlist.get(i).get(0) < indexlist.get(i + 1).get(0) && indexlist.get(i).get(0) > x.get(0)) {
                points.add(indexlist.get(i - 1));
            }
        }
        points.add(indexlist.get(indexlist.size() - 1));
        return points;
    }

    public static String typeAlongPath(List<IntTuple> path, String type) {
        String type_along_path = "";
        int index_type;
        for (int i = 1; i < path.size(); i++) {
            if (path.get(i - 1).get(1) < path.get(i).get(1)) {
                index_type = path.get(i - 1).get(1);
                type_along_path = type_along_path.concat(type.substring(index_type, index_type + 1));
            } else if (path.get(i - 1).get(1) > path.get(i).get(1)) {
                index_type = path.get(i).get(1);
                if (type.charAt(index_type) == 'f') {
                    type_along_path = type_along_path.concat("g");
                } else if (type.charAt(index_type) == 'g') {
                    type_along_path = type_along_path.concat("f");
                }
            } else if (path.get(i - 1).get(0) < path.get(i).get(0)) {
                type_along_path = type_along_path.concat("f");
            } else if (path.get(i - 1).get(0) > path.get(i).get(0)) {
                System.out.println("Die Pfadindizes stimmen nicht in die homogene Richung!!");
            } else {
                System.out.println("Die Pfadindizes stimmen nicht");
            }
        }
        return type_along_path;
    }

}

public class ReaderInputFile {
    public List<List<Point>> timeseries;

    public ReaderInputFile(String filename) throws FileNotFoundException {
        Scanner sc = new Scanner(new File(filename));
        int num_points_per_window = sc.nextInt();
        List<Point> list = new ArrayList<>();
        while (sc.hasNextLine()) {
            String data = sc.nextLine();
            // System.out.println("data:" + data);
            Scanner lineScanner = new Scanner(data);
            lineScanner.useLocale(Locale.US);
            List<Double> point_list = new ArrayList<>();
            // System.out.println("scanner has next: " + lineScanner.hasNext());
            while (lineScanner.hasNextDouble()) {
                Double tmp = lineScanner.nextDouble();
                point_list.add(tmp);
                // System.out.println("line scanner " + tmp);
            }
            if (point_list.size() > 0) {
                Point point = new Point(point_list);
                list.add(point);
            }
            lineScanner.close();
        }

        // teest:
        // for(int i=0; i<list.size(); i++){
        // System.out.println(list.get(i).getX());
        // }
        System.out.println(list.size());
        System.out.println(num_points_per_window);

        int num_windows = list.size() / num_points_per_window;
        List<List<Point>> out = new ArrayList<>();
        int total_iterator = 0;
        for (int i = 0; i < num_windows; i++) {
            List<Point> window = new ArrayList<>();
            for (int j = 0; j < num_points_per_window; j++) {
                window.add(list.get(total_iterator));
                // window.add(list.get(i*num_points_per_window+j));
                total_iterator++;
            }
            out.add(window);
        }
        sc.close();
        timeseries = out;
    }

}

public class SimplicialComplexTimeseries extends SimplicialComplex {
    private static Logger log = LoggerFactory.getLogger(SimplicialComplex.class);

    // public static List<List<Simplex>>
    // computeTimeseriesVRcomplex(List<List<Point>> timeseries, double
    // filtrationValue, int maxDimension, BinomialCoeffTable binom_coeff){
    // }

    public static List<List<Simplex>> computeTimeseriesVRcomplex(List<List<Point>> timeseries,
            List<DistanceMatrix> distanceMatrices, double filtrationValue, int maxDimension,
            BinomialCoeffTable binom_coeff) {
        List<List<Simplex>> simplices = new ArrayList<>();
        int total_index = 0;
        if (timeseries.size() == 1) {
            simplices.add(new ArrayList<>());
            for (int i = 0; i < timeseries.get(0).size(); i++) {
                simplices.get(0).add(new Simplex(i, 0));
            }
            DistanceMatrix distanceMatrix = distanceMatrices.get(0); // distance matrix von allen Punkten K_t \cup
                                                                     // K_{t+1}
            // Initialize neighbor sets
            Int2ObjectOpenHashMap<IntOpenHashSet> baseVertices = new Int2ObjectOpenHashMap<>();
            for (int i = 0; i < distanceMatrix.rows; i++) {
                baseVertices.put(i, new IntOpenHashSet());
            }
            int[] nonzeros_rows = distanceMatrix.getNonZeroRows();
            for (int i : nonzeros_rows) {
                IntOpenHashSet upperNeighbors = new IntOpenHashSet();
                List<Simplex> local_simplices = new ArrayList<>();
                int[] nonzero_cols = distanceMatrix.getNonZeroRowEntries(i);
                for (int j : nonzero_cols) {
                    if (j <= i) {
                        continue;
                    }
                    // Check if edge will appear in the filtration
                    if (distanceMatrix.get(i, j) <= filtrationValue) {
                        upperNeighbors.add(j);
                        local_simplices.add(new Simplex(binom_coeff.computeIndex(i, j), 1));
                    }
                }
                baseVertices.put(i, upperNeighbors);
                simplices.get(0).addAll(local_simplices);
            }
            int[] keys = baseVertices.keySet().toIntArray();
            for (int key : keys) {
                List<Simplex> local_simplices = new ArrayList<>();
                int[] vertices = new int[maxDimension + 1];
                vertices[0] = key;
                addCofaces(local_simplices, baseVertices.get(key), vertices, 0, baseVertices, binom_coeff);
                simplices.get(0).addAll(local_simplices);
            }
        }
        for (int t = 0; t < timeseries.size() - 1; t++) {
            int index_shift = total_index;
            // add vertices of step t and t+1
            simplices.add(new ArrayList<>());
            for (int i = 0; i < timeseries.get(t).size(); i++) {
                simplices.get(t).add(new Simplex(total_index, 0));
                total_index++;
            }
            int tmp_index = total_index;
            for (int i = 0; i < timeseries.get(t + 1).size(); i++) {
                simplices.get(t).add(new Simplex(tmp_index, 0));
                tmp_index++;
            }
            // ab hier kopiert:
            DistanceMatrix distanceMatrix = distanceMatrices.get(t); // distance matrix von allen Punkten K_t \cup
                                                                     // K_{t+1}
            // Initialize neighbor sets
            Int2ObjectOpenHashMap<IntOpenHashSet> baseVertices = new Int2ObjectOpenHashMap<>();
            for (int i = 0; i < distanceMatrix.rows; i++) {
                baseVertices.put(i + index_shift, new IntOpenHashSet());
            }
            int[] nonzeros_rows = distanceMatrix.getNonZeroRows();
            for (int i : nonzeros_rows) {
                IntOpenHashSet upperNeighbors = new IntOpenHashSet();
                List<Simplex> local_simplices = new ArrayList<>();
                int[] nonzero_cols = distanceMatrix.getNonZeroRowEntries(i);
                for (int j : nonzero_cols) {
                    if (j <= i) {
                        continue;
                    }
                    // Check if edge will appear in the filtration
                    if (distanceMatrix.get(i, j) <= filtrationValue) {
                        upperNeighbors.add(j + index_shift);
                        local_simplices.add(new Simplex(binom_coeff.computeIndex(i + index_shift, j + index_shift), 1));
                    }
                }
                baseVertices.put(i + index_shift, upperNeighbors);
                simplices.get(t).addAll(local_simplices);
            }
            int[] keys = baseVertices.keySet().toIntArray();
            for (int key : keys) {
                List<Simplex> local_simplices = new ArrayList<>();
                int[] vertices = new int[maxDimension + 1];
                vertices[0] = key;
                addCofaces(local_simplices, baseVertices.get(key), vertices, 0, baseVertices, binom_coeff);
                simplices.get(t).addAll(local_simplices);
            }
        }
        return simplices;
    }

    // TODO: vertexToTimeseries einmal ausrechnen und übergeben
    public static int belongsTo(Simplex simplex, List<List<Point>> timeseries, BinomialCoeffTable binomial_coeff) {
        List<Integer> indexListSimplices = new ArrayList<>();
        int indexTS = 0;
        indexListSimplices.add(indexTS);
        for (List<Point> window : timeseries) {
            indexTS += window.size();
            indexListSimplices.add(indexTS);
        }
        int[] vertices = Simplex.get_simplex_vertices(simplex.getIndex(), simplex.getDimension(),
                indexListSimplices.get(indexListSimplices.size() - 1), binomial_coeff);
        int first_vertex = vertices[0];
        int lower_bound_incl = 0, upper_bound_excl = 0;
        int time = -1;
        for (int i = 0; i < indexListSimplices.size(); i++) {
            if (indexListSimplices.get(i) > first_vertex) {
                lower_bound_incl = indexListSimplices.get(i - 1);
                upper_bound_excl = indexListSimplices.get(i);
                time = i - 1;
                break;
            }
        }
        for (int i = 1; i < vertices.length; i++) {
            if (vertices[i] < lower_bound_incl || vertices[i] >= upper_bound_excl) {
                if (vertices[i] < lower_bound_incl) {
                    time = time - 1;
                }
                return 2 * time + 1;
            }
        }
        return 2 * time;
    }

    // Kopiert von Topcat
    // candidates sind upperNeighbors von vertices[0]
    public static void addCofaces(List<Simplex> simplices, IntOpenHashSet candidates, int[] vertices, int index,
            Int2ObjectOpenHashMap<IntOpenHashSet> local_baseVertices, BinomialCoeffTable binomial_coeff) {
        if (index >= 2 && index < vertices.length) {
            List<Integer> t_vertices = new ArrayList<>(index + 1);
            for (int i = 0; i < index + 1; i++) {
                t_vertices.add(vertices[i]);
            }
            simplices.add(new Simplex(binomial_coeff.computeIndex(t_vertices), t_vertices.size() - 1));
        }
        if (index < vertices.length - 1) {
            IntIterator iterator = candidates.iterator();
            while (iterator.hasNext()) {
                vertices[index + 1] = iterator.nextInt();
                IntOpenHashSet lower_candidates = intersect(candidates, local_baseVertices.get(vertices[index + 1]));
                addCofaces(simplices, lower_candidates, vertices, index + 1, local_baseVertices, binomial_coeff);
            }
        }
    }

    // Kopiert von Topcat
    public static IntOpenHashSet intersect(IntOpenHashSet h1, IntOpenHashSet h2) {
        IntOpenHashSet intersection = new IntOpenHashSet();
        IntIterator iterator = h2.iterator();
        while (iterator.hasNext()) {
            int n = iterator.nextInt();
            if (h1.contains(n)) {
                intersection.add(n);
            }
        }
        return intersection;
    }

    // ToDo binomialCoeff nicht immer übergeben, sondern ein mal ausrechnen!!!
    // ToDo ab hier anpassen an eigenen Code, den Rest nochmal Korrektur lesen

    public static OwnSimplexStorageStructure computeSimplexStream(List<List<Point>> timeseries,
            List<DistanceMatrix> distanceMatrices, List<Double> filtrationValues, int maxDimension) {
        // Pick out the largest filtration value in each direction
        int num_vertices = 0;
        for (List<Point> window : timeseries) {
            num_vertices += window.size();
        }
        BinomialCoeffTable binomialCoeffTable = new BinomialCoeffTable(num_vertices, maxDimension);
        double maxFiltrationValue = filtrationValues.get(filtrationValues.size() - 1);
        IntTuple gridSize = IntTuple.zeros(2);
        gridSize.set(1, distanceMatrices.size() * 2 + 1);
        gridSize.set(0, filtrationValues.size() - 1); // warum -1 ?????????, bzw bei der zigzag Richtung auch -1 ? ?
        List<List<Double>> filtrationValuesList = new ArrayList<>();
        filtrationValuesList.add(filtrationValues);
        OwnSimplexStorageStructure storageStructure = new OwnSimplexStorageStructure(filtrationValuesList, gridSize,
                maxDimension, num_vertices);
        log.debug("Starting to compute simplicial complex...");
        // Compute the Vietoris-Rips complex for the maximal filtration value
        List<List<Simplex>> simplices = computeTimeseriesVRcomplex(timeseries, distanceMatrices, maxFiltrationValue,
                maxDimension, binomialCoeffTable);

        int num_simplices = 0;
        for (List<Simplex> lst : simplices) {
            num_simplices = num_simplices + lst.size();
        }
        log.debug("Finished computing simplicial complex. (Computed " + num_simplices + " number of simplices.)");
        log.debug("Starting to compute filtrationValues for each simplex...");
        // Compute the filtration indices for each simplex.

        // test simplex liste:
        // List<List<Simplex>> test = new ArrayList<>();

        // int time =0;
        for (List<Simplex> simplexList : simplices) {
            // time++;

            // List<Simplex> t = new ArrayList<>(); // test

            for (Simplex simplex : simplexList) {
                List<Integer> filtrationIndexes = calcFiltrationIndexes(simplex, timeseries, distanceMatrices,
                        filtrationValues, binomialCoeffTable);
                // if (filtrationIndexes != null) {
                if (!list_contains_simplex(simplex,
                        storageStructure.getSimplicesAt(simplex.getDimension(), new IntTuple(filtrationIndexes)))) {
                    storageStructure.addElement(simplex, new IntTuple(filtrationIndexes));
                    // t.add(simplex);
                }
                // } else {
                // log.error("Failed to add simplex " + simplex + " to simplex storage
                // structure");
                // }
            }
            // test.add(t);
        }

        // unnötig, nur zur Ausgabe da:
        /*
         * System.out.println("Test ");
         * for (int i=0; i<test.size();i++){
         * for (int j=0; j<test.get(i).size(); j++){
         * Simplex tmp = test.get(i).get(j);
         * int[] vertices =
         * Simplex.get_simplex_vertices(tmp.getIndex(),tmp.getDimension(),num_vertices,
         * binomialCoeffTable);
         * for (int vertex : vertices) {
         * System.out.print(" " + vertex);
         * }
         * System.out.println();
         * }
         * }
         */

        log.debug("Finished computing filtrationValues.");
        return storageStructure;
    }

    public static List<Integer> calcFiltrationIndexes(Simplex simplex, List<List<Point>> timeseries,
            List<DistanceMatrix> distanceMatrices, List<Double> filtrationValues, BinomialCoeffTable binomial_coeff) {
        // Find maximum value of the weights on the edges for each metric
        int num_vertices = 0;
        for (List<Point> window : timeseries) {
            num_vertices += window.size();
        }
        int[] vertices = Simplex.get_simplex_vertices(simplex.getIndex(), simplex.getDimension(), num_vertices,
                binomial_coeff);
        int time = belongsTo(simplex, timeseries, binomial_coeff);

        List<Integer> filtrationIndices = new ArrayList<>();
        filtrationIndices.add(time); // time stimmt wahrscheinlich nicht ganz!!!!!!!!!!!!!!! time ist der index des
                                     // zigzags -> soll ja auch so sein
        int index_shift;
        DistanceMatrix distanceMatrix;
        if (time % 2 == 0 && time / 2 < distanceMatrices.size()) {
            distanceMatrix = distanceMatrices.get(time / 2);
            index_shift = time / 2 * timeseries.get(0).size();
        } else if (time % 2 == 0) {
            distanceMatrix = distanceMatrices.get(time / 2 - 1);
            index_shift = (time / 2 - 1) * timeseries.get(0).size();
        } else {
            distanceMatrix = distanceMatrices.get(time / 2); // time/2+1
            index_shift = (time / 2) * timeseries.get(0).size();
        }
        double f_max = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < vertices.length; i++) {
            for (int j = i; j < vertices.length; j++) {
                double f = distanceMatrix.get(vertices[i] - index_shift, vertices[j] - index_shift);
                if (f > f_max) {
                    f_max = f;
                }
            }
        }
        int filtrationIndex = 0, i = 0;
        while (i < filtrationValues.size() && f_max > filtrationValues.get(i)) {
            filtrationIndex = ++i;
        }
        filtrationIndices.add(filtrationIndex);
        return new ArrayList<>(Arrays.asList(filtrationIndices.get(1), filtrationIndices.get(0))); // reversion of
                                                                                                   // indices
    }

    public static boolean list_contains_simplex(Simplex simplex, List<Simplex> list) {
        if (list == null) {
            return false;
        }
        for (Simplex simpl_comp : list) {
            if (simplex.equals(simpl_comp)) {
                return true;
            }
        }
        return false;
    }
}

public class EntryPoint {
    private Filtration filtration = null;
    private OwnSimplexStorageStructure ownSimplexStorageStructure = null;
    private Path path = null;
    private ReaderInputFile readerInputFile = null;
    private SimplicialComplexTimeseries simplicialComplexTimeseries = null;

    public EntryPoint() {
    }

    public Filtration getClassFiltration(List<List<Point>> timeseries, List<DistanceMatrix> distanceMatrices,
            List<Double> filtrationValues, int maxDimension) {
        if (filtration == null) {
            filtration = new Filtration(timeseries, distanceMatrices, filtrationValues, maxDimension);
        }
        return filtration;
    }

    public OwnSimplexStorageStructure getClassOwnSimplexStorageStructure(List<List<Point>> timeseries,
            List<DistanceMatrix> distanceMatrices, List<Double> filtrationValues, int maxDimension) {
        if (ownSimplexStorageStructure == null) {
            ownSimplexStorageStructure = SimplicialComplexTimeseries.computeSimplexStream(timeseries, distanceMatrices,
                    filtrationValues, maxDimension);
        }
        return ownSimplexStorageStructure;
    }

    public Path getClassPath() {
        if (path == null) {
            path = new Path();
        }
        return path;
    }

    public ReaderInputFile getClassReaderInputFile(String filename) {
        if (readerInputFile == null) {
            try {
                readerInputFile = new ReaderInputFile(filename);
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
        }
        return readerInputFile;
    }

    public SimplicialComplexTimeseries getClassSimplicialComplexTimeseries() {
        if (simplicialComplexTimeseries == null) {
            simplicialComplexTimeseries = new SimplicialComplexTimeseries();
        }
        return simplicialComplexTimeseries;
    }

    public static void main(String[] args) {
        EntryPoint entryPoint = new EntryPoint();
        GatewayServer server = new GatewayServer(entryPoint);
        server.start();
        System.out.println("Gateway Server Started");
    }
}