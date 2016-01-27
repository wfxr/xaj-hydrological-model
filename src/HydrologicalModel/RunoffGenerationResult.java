package HydrologicalModel;

/**
 * 产流计算结果
 *
 * Created by Wenxuan on 2016/1/11.
 */
public class RunoffGenerationResult {
    /**
     * 序列长度
     */
    public final int Length;

    /**
     * 序列长度
     */
    public double[] PE;

    /**
     * 产流量序列
     */
    public double[] R;

    /**
     * 土壤含水量序列
     */
    public double[] W;

    /**
     * @param n 序列长度
     */
    public RunoffGenerationResult(int n) {
        PE = new double[n];
        R = new double[n];
        W = new double[n];
        Length = n;
    }

    /**
     * 设置序列指定位置处的产流计算结果
     *
     * @param idx 序列位置索引
     * @param pe 净雨量
     * @param r 产流量
     * @param w 土壤含水量
     */
    public void Set(int idx, double pe, double r, double w){
        PE[idx] = pe;
        R[idx] = r;
        W[idx] = w;
    }
}
