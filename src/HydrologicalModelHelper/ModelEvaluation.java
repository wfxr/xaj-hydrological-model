package HydrologicalModelHelper;

import java.util.stream.DoubleStream;

import static HydrologicalModelHelper.ParameterValidation.*;

/**
 * 模型评价工具类
 *
 * Created by Wenxuan on 1/27/2016.
 * Email: wenxuan-zhang@outlook.com
 */
public class ModelEvaluation {
    /**
     * 计算模型预测结果的纳什效率系数
     *
     * @param obs 实际观测值序列
     * @param pre 模型预测值序列
     * @return 纳什效率系数
     */
    public static double NashSutcliffeEfficiency(double[] obs, double[] pre) {
        int n = CheckIdenticalLength(obs, pre);
        double avg = DoubleStream.of(obs).average().getAsDouble();
        double s1 = 0, s2 = 0;
        for (int i = 0; i < n; ++i) {
            s1 += Math.pow(obs[i] - pre[i], 2);
            s2 += Math.pow(obs[i] - avg, 2);
        }
        return 1 - s1 / s2;
    }
}
