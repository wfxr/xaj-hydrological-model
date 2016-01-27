package HydrologicalModelHelper;

import static com.google.common.base.Preconditions.*;

/**
 * Created by Wenxuan on 1/27/2016.
 * Email: wenxuan-zhang@outlook.com
 */
public class ParameterValidation {
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
