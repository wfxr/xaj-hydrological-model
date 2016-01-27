package HydrologicalModelHelper;

import static com.google.common.base.Preconditions.*;

/**
 * 参数验证工具类
 *
 * Created by Wenxuan on 1/27/2016.
 * Email: wenxuan-zhang@outlook.com
 */
public class ParameterValidation {
    /**
     * 确保所有的待检验的数组长度一致
     *
     * @param arrays 待检验的数组
     * @return 如果所有数组长度一致则返回数组的长度
     */
    public static int CheckIdenticalLength(double[]... arrays) {
        checkNotNull(arrays);
        checkElementIndex(0, arrays.length);
        boolean identical = true;
        for (int i = 1; i < arrays.length; ++i)
            if (arrays[i].length != arrays[0].length)
                identical = false;
        checkArgument(identical, "序列长度不一致");
        return arrays[0].length;
    }
}
