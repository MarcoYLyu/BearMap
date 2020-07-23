/**
 * This class provides all code necessary to take a query box and produce
 * a query result. The getMapRaster method must return a Map containing all
 * seven of the required fields, otherwise the front end code will probably
 * not draw the output correctly.
 */
public class Rasterer {
    /**
     * The max image depth level.
     */
    public static final int MAX_DEPTH = 7;
    public static final String IMGTEMPLATE = "d%d_x%d_y%d.png";

    public static void main(String[] args) {
        Rasterer ras = new Rasterer();

    }

    /**
     * Takes a user query and finds the grid of images that best matches the query. These images
     * will be combined into one big image (rastered) by the front end. The grid of images must obey
     * the following properties, where image in the grid is referred to as a "tile".
     * <ul>
     * <li>The tiles collected must cover the most longitudinal distance per pixel (LonDPP)
     * possible, while still covering less than or equal to the amount of longitudinal distance
     * per pixel in the query box for the user viewport size.</li>
     * <li>Contains all tiles that intersect the query bounding box that fulfill the above
     * condition.</li>
     * <li>The tiles must be arranged in-order to reconstruct the full image.</li>
     * </ul>
     *
     * @param params The RasterRequestParams containing coordinates of the query box and the browser
     *               viewport width and height.
     * @return A valid RasterResultParams containing the computed results.
     */
    public RasterResultParams getMapRaster(RasterRequestParams params) {
        if (params.lrlat >= params.ullat || params.ullon > params.lrlon) {
            return RasterResultParams.queryFailed();
        }

        int depth = getDepth(params.lrlon, params.ullon, params.w);

        int ulx = getX(params.ullon, depth);
        int uly = getY(params.ullat, depth);
        int lrx = getX(params.lrlon, depth);
        int lry = getY(params.lrlat, depth);

        String[][] renderGrid = getRenderGrid(ulx, uly, lrx, lry, depth);
        double rasterUlLon = MapServer.ROOT_ULLON + getUnitLon(depth) * ulx;
        double rasterUlLat = MapServer.ROOT_ULLAT - getUnitLat(depth) * uly;
        double rasterLrLon = MapServer.ROOT_ULLON + getUnitLon(depth) * (lrx + 1);
        double rasterLrLat = MapServer.ROOT_ULLAT - getUnitLat(depth) * (lry + 1);

        RasterResultParams.Builder builder = new RasterResultParams.Builder();
        builder.setRenderGrid(renderGrid);
        builder.setRasterLrLat(rasterLrLat);
        builder.setRasterLrLon(rasterLrLon);
        builder.setRasterUlLat(rasterUlLat);
        builder.setRasterUlLon(rasterUlLon);
        builder.setDepth(depth);
        builder.setQuerySuccess(true);

        return builder.create();
    }

    private double getDepthDPP(int n) {
        return lonDPP(MapServer.ROOT_LRLON, MapServer.ROOT_ULLON, 256 * Math.pow(2, n));
    }

    /* Private method */
    public int getDepth(double lrlon, double ullon, double width) {
        double target = lonDPP(lrlon, ullon, width);
        for (int i = 0; i <= MAX_DEPTH; ++i) {
            if (getDepthDPP(i) <= target) {
                return i;
            }
        }
        return MAX_DEPTH;
    }

    /* Private method */
    private double getUnitLon(int n) {
        return (MapServer.ROOT_LRLON - MapServer.ROOT_ULLON) / Math.pow(2, n);
    }

    private double getUnitLat(int n) {
        return (MapServer.ROOT_ULLAT - MapServer.ROOT_LRLAT) / Math.pow(2, n);
    }

    public int getX(double x, int n) {
        if (x < MapServer.ROOT_ULLON) {
            return 0;
        } else if (x > MapServer.ROOT_LRLON) {
            return (int) (Math.pow(2, n)) - 1;
        }
        return (int) ((x - MapServer.ROOT_ULLON) / getUnitLon(n));
    }

    public int getY(double y, int n) {
        if (y > MapServer.ROOT_ULLAT) {
            return 0;
        } else if (y < MapServer.ROOT_LRLAT) {
            return (int) (Math.pow(2, n)) - 1;
        }
        return (int) ((MapServer.ROOT_ULLAT - y) / getUnitLat(n));
    }

    public String[][] getRenderGrid(int ulx, int uly, int lrx, int lry, int n) {
        String[][] renderGrid = new String[lry - uly + 1][lrx - ulx + 1];
        for (int i = 0; i < lry - uly + 1; i++) {
            for (int j = 0; j < lrx - ulx + 1; j++) {
                renderGrid[i][j] = String.format(IMGTEMPLATE, n, ulx + j, uly + i);
            }
        }
        return renderGrid;
    }

    /**
     * Calculates the lonDPP of an image or query box
     *
     * @param lrlon Lower right longitudinal value of the image or query box
     * @param ullon Upper left longitudinal value of the image or query box
     * @param width Width of the query box or image
     * @return lonDPP
     */
    private double lonDPP(double lrlon, double ullon, double width) {
        return (lrlon - ullon) / width;
    }
}
