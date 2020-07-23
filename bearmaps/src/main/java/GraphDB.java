import org.xml.sax.SAXException;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * Graph for storing all of the intersection (vertex) and road (edge) information.
 * Uses your GraphBuildingHandler to convert the XML files into a graph. Your
 * code must include the vertices, adjacent, distance, closest, lat, and lon
 * methods. You'll also need to include instance variables and methods for
 * modifying the graph (e.g. addNode and addEdge).
 *
 * @author Kevin Lowe, Antares Chen, Kevin Lin
 */
public class GraphDB {
    /**
     * Radius of the Earth in miles.
     */
    private static final int R = 3963;
    /**
     * Latitude centered on Berkeley.
     */
    private static final double ROOT_LAT = (MapServer.ROOT_ULLAT + MapServer.ROOT_LRLAT) / 2;
    /**
     * Longitude centered on Berkeley.
     */
    private static final double ROOT_LON = (MapServer.ROOT_ULLON + MapServer.ROOT_LRLON) / 2;
    /**
     * Scale factor at the natural origin, Berkeley. Prefer to use 1 instead of 0.9996 as in UTM.
     *
     * @source https://gis.stackexchange.com/a/7298
     */
    private static final double K0 = 1.0;
    /**
     * idToNode: map from id to node/vertex
     * idToEdge: map from wayID to edge
     */
    Map<Long, Node> idToNode;
    Map<Long, Edge> idToEdge;
    KDTree tree;

    /**
     * This constructor creates and starts an XML parser, cleans the nodes, and prepares the
     * data structures for processing. Modify this constructor to initialize your data structures.
     *
     * @param dbPath Path to the XML file to be parsed.
     */
    public GraphDB(String dbPath) {
        idToNode = new HashMap<>();
        idToEdge = new HashMap<>();

        File inputFile = new File(dbPath);
        try (FileInputStream inputStream = new FileInputStream(inputFile)) {
            SAXParserFactory factory = SAXParserFactory.newInstance();
            SAXParser saxParser = factory.newSAXParser();
            saxParser.parse(inputStream, new GraphBuildingHandler(this));
        } catch (ParserConfigurationException | SAXException | IOException e) {
            e.printStackTrace();
        }
        clean();
        tree = getKDTree();
    }

    /**
     * Helper to process strings into their "cleaned" form, ignoring punctuation and capitalization.
     *
     * @param s Input string.
     * @return Cleaned string.
     */
    private static String cleanString(String s) {
        return s.replaceAll("[^a-zA-Z ]", "").toLowerCase();
    }

    /**
     * Return the Euclidean x-value for some point, p, in Berkeley. Found by computing the
     * Transverse Mercator projection centered at Berkeley.
     *
     * @param lon The longitude for p.
     * @param lat The latitude for p.
     * @return The flattened, Euclidean x-value for p.
     * @source https://en.wikipedia.org/wiki/Transverse_Mercator_projection
     */
    static double projectToX(double lon, double lat) {
        double dlon = Math.toRadians(lon - ROOT_LON);
        double phi = Math.toRadians(lat);
        double b = Math.sin(dlon) * Math.cos(phi);
        return (K0 / 2) * Math.log((1 + b) / (1 - b));
    }

    /**
     * Return the Euclidean y-value for some point, p, in Berkeley. Found by computing the
     * Transverse Mercator projection centered at Berkeley.
     *
     * @param lon The longitude for p.
     * @param lat The latitude for p.
     * @return The flattened, Euclidean y-value for p.
     * @source https://en.wikipedia.org/wiki/Transverse_Mercator_projection
     */
    static double projectToY(double lon, double lat) {
        double dlon = Math.toRadians(lon - ROOT_LON);
        double phi = Math.toRadians(lat);
        double con = Math.atan(Math.tan(phi) / Math.cos(dlon));
        return K0 * (con - Math.toRadians(ROOT_LAT));
    }

    public KDTree getKDTree() {
        if (tree == null) {
            List<Long> ids = new ArrayList<>(idToNode.keySet());
            tree = new KDTree(ids);
//            for (long id : idToNode.keySet()) {
//                tree.root = tree.insertKD(tree.root, id, 0);
//            }
        }
        return tree;
    }

    /**
     * @param id  id of added Node
     * @param lat lat of added Node
     * @param lon lon og added Node
     */
    public void addNode(long id, double lat, double lon) {
        idToNode.put(id, new Node(id, lat, lon));
    }

    /**
     * @param wayID          id of way
     * @param posConnections all nodes on edge
     */
    public void addEdge(long wayID, LinkedList<Long> posConnections) {
        idToEdge.put(wayID, new Edge(wayID, posConnections));

        for (int i = 0; i < posConnections.size() - 1; i++) {
            long idToConnect1 = posConnections.get(i);
            long idToConnect2 = posConnections.get(i + 1);
            Node n1 = idToNode.get(idToConnect1);
            Node n2 = idToNode.get(idToConnect2);
            n1.becomeNeighbors(n2);
        }
    }

    public void removeNode(long id) {
        idToNode.remove(id);
    }

    /**
     * Remove nodes with no connections from the graph.
     * While this does not guarantee that any two nodes in the remaining graph are connected,
     * we can reasonably assume this since typically roads are connected.
     */
    private void clean() {
        HashSet<Long> idsToRemove = new HashSet<>();
        for (Node n : idToNode.values()) {
            if (n.adjList.isEmpty()) {
                idsToRemove.add(n.id);
            }
        }
        for (long id : idsToRemove) {
            idToNode.remove(id);
        }
    }

    /**
     * Returns the longitude of vertex <code>v</code>.
     *
     * @param v The ID of a vertex in the graph.
     * @return The longitude of that vertex, or 0.0 if the vertex is not in the graph.
     */
    double lon(long v) {
        if (idToNode.containsKey(v)) {
            Node n = idToNode.get(v);
            return n.lon;
        }
        return 0;
    }

    /**
     * Returns the latitude of vertex <code>v</code>.
     *
     * @param v The ID of a vertex in the graph.
     * @return The latitude of that vertex, or 0.0 if the vertex is not in the graph.
     */
    double lat(long v) {
        if (idToNode.containsKey(v)) {
            Node n = idToNode.get(v);
            return n.lat;
        }
        return 0;
    }

    /**
     * Returns the x-axis of vertex <code>v</code>.
     *
     * @param v The ID of a vertex in the graph.
     * @return The x of that vertex, or 0.0 if the vertex is not in the graph.
     */
    double getX(long v) {
        if (idToNode.containsKey(v)) {
            Node n = idToNode.get(v);
            return n.x;
        }
        return 0;
    }

    /**
     * Returns the Y-axis of vertex <code>v</code>.
     *
     * @param v The ID of a vertex in the graph.
     * @return The x of that vertex, or 0.0 if the vertex is not in the graph.
     */
    double getY(long v) {
        if (idToNode.containsKey(v)) {
            Node n = idToNode.get(v);
            return n.y;
        }
        return 0;
    }
    /**
     * Returns an iterable of all vertex IDs in the graph.
     *
     * @return An iterable of all vertex IDs in the graph.
     */
    Iterable<Long> vertices() {
        return idToNode.keySet();
    }

    /**
     * Returns an iterable over the IDs of all vertices adjacent to <code>v</code>.
     *
     * @param v The ID for any vertex in the graph.
     * @return An iterable over the IDs of all vertices adjacent to <code>v</code>, or an empty
     * iterable if the vertex is not in the graph.
     */
    Iterable<Long> adjacent(long v) {
        if (idToNode.containsKey(v)) {
            Node n = idToNode.get(v);
            return n.adjList;
        }
        return Collections.emptySet();
    }

    /**
     * Returns the great-circle distance between two vertices, v and w, in miles.
     * Assumes the lon/lat methods are implemented properly.
     *
     * @param v The ID for the first vertex.
     * @param w The ID for the second vertex.
     * @return The great-circle distance between vertices and w.
     * @source https://www.movable-type.co.uk/scripts/latlong.html
     */
    public double distance(long v, long w) {
        double phi1 = Math.toRadians(lat(v));
        double phi2 = Math.toRadians(lat(w));
        double dphi = Math.toRadians(lat(w) - lat(v));
        double dlambda = Math.toRadians(lon(w) - lon(v));

        double a = Math.sin(dphi / 2.0) * Math.sin(dphi / 2.0);
        a += Math.cos(phi1) * Math.cos(phi2) * Math.sin(dlambda / 2.0) * Math.sin(dlambda / 2.0);
        double c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
        return R * c;
    }

    /**
     * Returns the ID of the vertex closest to the given longitude and latitude.
     *
     * @param lon The given longitude.
     * @param lat The given latitude.
     * @return The ID for the vertex closest to the <code>lon</code> and <code>lat</code>.
     */
    public long closest2(double lon, double lat) {
        addNode(-1, lat, lon);
        long closestId = 0;
        double closestDist = Double.MAX_VALUE;

        for (long n : idToNode.keySet()) {
            double dist = distance(n, -1);
            if (dist < closestDist && n != -1) {
                closestId = n;
                closestDist = dist;
            }
        }
        removeNode(-1);
        return closestId;
    }

    public long closest(double lon, double lat) {
        addNode(-1, lat, lon);
//        long id = tree.searchClosest(null, null, tree.root, 0, Double.MAX_VALUE).nodeID;
        double[] coord = tree.getXYCoord(-1);
        double[] bestCoord = tree.getXYCoord(tree.root.nodeID);

        tree.setminDistanceToTarget(euclidean2(bestCoord[0], coord[0], bestCoord[1], coord[1]));
        long id = tree.searchClosest2(tree.root, null, tree.root, 0, coord).nodeID;
        removeNode(-1);
        return id;
    }

    /**
     * In linear time, collect all the names of OSM locations that prefix-match the query string.
     *
     * @param prefix Prefix string to be searched for. Could be any case, with our without
     *               punctuation.
     * @return A <code>List</code> of the full names of locations whose cleaned name matches the
     * cleaned <code>prefix</code>.
     */
    public List<String> getLocationsByPrefix(String prefix) {
        return Collections.emptyList();
    }

    /**
     * Collect all locations that match a cleaned <code>locationName</code>, and return
     * information about each node that matches.
     *
     * @param locationName A full name of a location searched for.
     * @return A <code>List</code> of <code>LocationParams</code> whose cleaned name matches the
     * cleaned <code>locationName</code>
     */
    public List<LocationParams> getLocations(String locationName) {
        return Collections.emptyList();
    }

    /**
     * Returns the initial bearing between vertices <code>v</code> and <code>w</code> in degrees.
     * The initial bearing is the angle that, if followed in a straight line along a great-circle
     * arc from the starting point, would take you to the end point.
     * Assumes the lon/lat methods are implemented properly.
     *
     * @param v The ID for the first vertex.
     * @param w The ID for the second vertex.
     * @return The bearing between <code>v</code> and <code>w</code> in degrees.
     * @source https://www.movable-type.co.uk/scripts/latlong.html
     */
    double bearing(long v, long w) {
        double phi1 = Math.toRadians(lat(v));
        double phi2 = Math.toRadians(lat(w));
        double lambda1 = Math.toRadians(lon(v));
        double lambda2 = Math.toRadians(lon(w));

        double y = Math.sin(lambda2 - lambda1) * Math.cos(phi2);
        double x = Math.cos(phi1) * Math.sin(phi2);
        x -= Math.sin(phi1) * Math.cos(phi2) * Math.cos(lambda2 - lambda1);
        return Math.toDegrees(Math.atan2(y, x));
    }

    static double euclidean2(double x1, double x2, double y1, double y2) {
        return Math.pow(x1 - x2, 2) + Math.pow(y1 - y2, 2);
    }

    private class Node {
        LinkedList<Long> adjList;
        private long id;
        private double lon;
        private double lat;

        private double x;
        private double y;

        Node(long id, double lat, double lon) {
            this.id = id;
            this.lon = lon;
            this.lat = lat;
            this.x = projectToX(lon, lat);
            this.y = projectToY(lon, lat);
            adjList = new LinkedList<>();
        }

        private void becomeNeighbors(Node n) {
            adjList.add(n.id);
            n.adjList.add(id);
        }
    }

    private class Edge {
        private long id;
        private LinkedList<Long> posConnections;

        Edge(long id, LinkedList<Long> posConnections) {
            this.id = id;
            this.posConnections = posConnections;
        }
    }

    private class KDTree {

        private KDNode root;

        private double minDistanceToTarget;

        KDTree() {

        }

        KDTree(List<Long> nodeIds) {
            root = buildKDTree(nodeIds, 0);
        }

        public void setminDistanceToTarget(double dMin) {
            this.minDistanceToTarget = dMin;
        }


        public KDNode buildKDTree(List<Long> nodeIds, int depth) {
//            Iterator<Long> itr = nodeIds.iterator();
//            while (itr.hasNext()) {
//                root = insertKD(root, itr.next(), 0);
//            }
            if (nodeIds.isEmpty()) {
                return null;
            }
            int cd = depth % 2;
            if (cd == 0) {
                nodeIds.sort((x, y) ->
                        new Double(projectToX(lon(x), lat(x)))
                                .compareTo(projectToX(lon(y), lat(y))));
            } else {
                nodeIds.sort((x, y) ->
                        new Double(projectToY(lon(x), lat(x)))
                                .compareTo(projectToY(lon(y), lat(y))));
            }

            int median = nodeIds.size() / 2;

            return new KDNode(nodeIds.get(median),
                    buildKDTree(nodeIds.
                            subList(0, median), depth + 1),
                    buildKDTree(nodeIds.
                            subList(median + 1, nodeIds.size()), depth + 1));
        }

        private double[] getXYCoord(long nodeID) {
            double[] result = {getX(nodeID), getY(nodeID)};
            return result;
        }

        public KDNode searchClosest2(KDNode best, KDNode pre, KDNode n, int depth, double[] coord) {
            if (n == null) {
                return best;
            }
            double[] curCoord = getXYCoord(n.nodeID);
            double curToTarget = euclidean2(curCoord[0], coord[0], curCoord[1], coord[1]);
            if (curToTarget < this.minDistanceToTarget) {
                this.setminDistanceToTarget(curToTarget);
                best = n;
            }
            int cd = depth % 2;
            double checkToTarget = Math.pow(curCoord[cd] - coord[cd], 2);

            if (coord[cd] < curCoord[cd]) {  //check left first
                best = searchClosest2(best, n, n.left, depth + 1, coord);
                if (this.minDistanceToTarget > checkToTarget) {
                    best = searchClosest2(best, n, n.right, depth + 1, coord);
                }
            } else {
                best = searchClosest2(best, n, n.right, depth + 1, coord);
                if (this.minDistanceToTarget > checkToTarget) {
                    best = searchClosest2(best, n, n.left, depth + 1, coord);
                }
            }
            return best;
        }

        private class KDNode {
            long nodeID;
            KDNode left;
            KDNode right;

            KDNode(long id, KDNode l, KDNode r) {
                nodeID = id;
                left = l;
                right = r;
            }
        }
    }
}
