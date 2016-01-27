package HydrologicalModel;

/**
 * 水源划分计算结果
 *
 * Created by Wenxuan on 2016/1/11.
 */
@SuppressWarnings("SpellCheckingInspection")
public class SourcePartitionResult {
    /**
     * 序列长度
     */
    public final int Length;

    /**
     * 不透水面积产流序列
     */
    public double[] RIMP;

    /**
     * 地表径流产流序列
     */
    public double[] RS;

    /**
     * 壤中流产流序列
     */
    public double[] RSS;

    /**
     * 地下径流产流序列
     */
    public double[] RG;

    /**
     * 流域平均自由含水量序列
     */
    public double[] S;

    /**
     * @param n 序列长度
     */
    public SourcePartitionResult(int n) {
        RIMP = new double[n];
        RS = new double[n];
        RSS = new double[n];
        RG = new double[n];
        S = new double[n];
        Length = n;
    }

    /**
     * 设置序列指定位置处的水源划分计算结果
     *
     * @param idx 序列位置索引
     * @param rimp 不透水面积产流
     * @param rs 地表径流产流
     * @param rss 壤中流产流
     * @param rg 地下径流产流
     * @param s 流域平均自由含水量
     */
    public void Set(int idx, double rimp, double rs, double rss, double rg, double s){
        RIMP[idx] = rimp;
        RS[idx] = rs;
        RSS[idx] = rss;
        RG[idx] = rg;
        S[idx] = s;
    }
}
